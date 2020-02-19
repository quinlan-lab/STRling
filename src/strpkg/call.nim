import kmer
import math
import algorithm
import msgpack4nim
import strutils
import times
import random
import tables
import hts/bam
import hts/fai
import strformat
import math
import argparse

import ./cluster
import ./collect
import ./utils
import ./genotyper
import ./genome_strs
import ./extract
import ./callclusters
import ./unpack

export tread
export Soft

const n_sims = 100
const each_n = 2
proc string_hash(s:string, sim_i:int): uint64 {.inline.} =
  # djb's string hash
  result = 5381
  for c in s:
    result = ((result shl 5) + result) + c.uint64

  for c in $sim_i:
    result = ((result shl 5) + result) + c.uint64

proc filter_sims[T:tread|Support](trs: seq[T], sim_i: int, each_n:int): seq[T] =
  result = trs
  if result.len == 0: return
  randomize(sim_i + 1)
  shuffle(result)
  result = result[0..<max(1, int(result.len/each_n))]

proc call_main*() =
  var p = newParser("strling call"):
    option("-f", "--fasta", help="path to fasta file")
    option("-m", "--min-support", help="minimum number of supporting reads for a locus to be reported", default="6")
    option("-c", "--min-clip", help="minimum number of supporting clipped reads for each side of a locus", default="0")
    option("-t", "--min-clip-total", help="minimum total number of supporting clipped reads for a locus", default="0")
    option("-q", "--min-mapq", help="minimum mapping quality (does not apply to STR reads)", default="40")
    option("-l", "--loci", help="Annoated bed file specifying additional STR loci to genotype. Format is: chr start stop repeatunit [name]")
    option("-b", "--bounds", help="STRling -bounds.txt file (usually produced by strling merge) specifying additional STR loci to genotype.")
    option("-o", "--output-prefix", help="prefix for output files", default="strling")
    flag("-v", "--verbose")
    arg("bam", help="path to bam file")
    arg("bin", help="bin file previously created by `strling extract`")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "call":
    argv = argv[1..argv.high]

  if len(argv) == 0: argv = @["-h"]

  var args = p.parse(argv)
  if args.help:
    quit 0

  var t0 = cpuTime()
  var ibam_dist:Bam
  var min_support = parseInt(args.min_support)
  var min_clip = uint16(parseInt(args.min_clip))
  var min_clip_total = uint16(parseInt(args.min_clip_total))
  var min_mapq = uint8(parseInt(args.min_mapq))
  var skip_reads = 100000

  if args.loci != "":
    if not fileExists(args.loci):
      quit "couldn't open loci file"

  if args.bounds != "":
    if not fileExists(args.bounds):
      quit "couldn't open bounds file"

  if not open(ibam_dist, args.bam, fai=args.fasta, threads=0):
    quit "couldn't open bam"

  var cram_opts = 8191 - SAM_RNAME.int - SAM_RGAUX.int - SAM_QUAL.int - SAM_SEQ.int
  discard ibam_dist.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, cram_opts)

  var frag_dist = ibam_dist.fragment_length_distribution(skip_reads=skip_reads)
  #echo frag_dist[0..<1000]
  var frag_median = frag_dist.median
  if args.verbose:
    stderr.write_line "Calculated median fragment length:", frag_median
    stderr.write_line "10th, 90th percentile of fragment length:", $frag_dist.median(0.1), " ", $frag_dist.median(0.9)

  ibam_dist.close()
  ibam_dist = nil
  var ibam:Bam
  if not open(ibam, args.bam, fai=args.fasta, threads=0, index=true):
    quit "couldn't open bam"
  cram_opts = 8191 - SAM_RGAUX.int - SAM_QUAL.int
  discard ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, cram_opts)

  var opts = Options(median_fragment_length: frag_median,
                      min_clip: min_clip, min_clip_total: min_clip_total,
                      min_support: min_support, min_mapq: min_mapq,
                      window: frag_dist.median(0.98),
                      targets: ibam.hdr.targets)

  # Unpack STR reads from bin file and put them in a table by repeat unit and chromosome
  var treads_by_tid_rep = newTable[tid_rep, seq[tread]](8192)

  var fs = newFileStream(args.bin, fmRead)
  var extracted = fs.unpack_file()
  doAssert extracted.targets.same(ibam.hdr.targets)
  fs.close()
  for r in extracted.reads:
    treads_by_tid_rep.mgetOrPut((r.tid, r.repeat), newSeq[tread]()).add(r)
  # TODO: Check all bin file came from the current version of STRling
  for k, trs in treads_by_tid_rep.mpairs:
    trs.sort(proc(a, b:tread): int =
      cmp(a.position, b.position)
    )

  ### discovery
  var
    gt_fh:File
    bounds_fh:File
    unplaced_fh:File
  if not open(gt_fh, args.output_prefix & "-genotype.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(bounds_fh, args.output_prefix & "-bounds.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(unplaced_fh, args.output_prefix & "-unplaced.txt", mode=fmWrite):
    quit "couldn't open output file"

  # Write headers
  bounds_fh.write_line(bounds_header & "\tdepth")
  gt_fh.write_line(gt_header)

  when defined(debug):
    var
      reads_fh:File
      span_fh:File
    if not open(reads_fh, args.output_prefix & "-reads.txt", mode=fmWrite):
      quit "couldn't open output file"
    if not open(span_fh, args.output_prefix & "-spanning.txt", mode=fmWrite):
      quit "couldn't open output file"
    reads_fh.write_line &"#chrom\tpos\tstr\tsoft_clip\tstr_count\tqname\tcluster_id"

  var loci: seq[Bounds]
  if args.loci != "":
    # Parse bed file of regions. These will also be genotyped
    loci = parse_bed(args.loci, opts.targets, opts.window.uint32)
    stderr.write_line &"Read {len(loci)} loci from {args.loci}"

  var bounds: seq[Bounds]
  if args.bounds != "":
    # Parse -bounds.txt file of regions. These will also be genotyped
    bounds = parse_bounds(args.bounds, opts.targets)
    stderr.write_line &"Read {len(bounds)} bounds from {args.bounds}"
  else:
    quit "-b bounds file is required for this branch"

  # Merge loci and bounds, with loci overwriting bounds
  for bound in bounds.mitems:
    for i, locus in loci:
      if locus.overlaps(bound):
        bound.name = locus.name
        bound.left = locus.left
        bound.right = locus.right
        # Remove locus from loci (therefore will use first matching bound and locus)
        loci.del(i)
        break

  # Add any remaining loci to bounds
  for locus in loci:
    bounds.add(locus)


  # Assign STR reads to provided bounds and genotype them
  t0 = cpuTime()

  var spansByBound = initTable[Bounds, seq[Support]]()
  var treadsByBound = initTable[Bounds, seq[tread]]()
  var medByBound = initTable[Bounds, int]()
  for sim_i in 0..<n_sims:
    echo "sim:", sim_i, " time:", cpuTime() - t0
    t0 = cpuTime()
    var unplaced_counts = initCountTable[string]()
    var genotypes_by_repeat = initTable[string, seq[Call]]()

    for bound in bounds.mitems:

      if bound.right - bound.left > 1000'u32:
        stderr.write_line "large bounds:" & $bound & " skipping"
        continue

      if bound notin treadsByBound:
        treadsByBound[bound] = assign_reads_locus(bound, treads_by_tid_rep)
      var str_reads = treadsByBound[bound].filter_sims(sim_i, each_n)

      if bound notin spansByBound:
        var (s, md) = ibam.spanners(bound, opts.window, frag_dist, opts.min_mapq)
        spansByBound[bound] = s
        medByBound[bound] = md

      var spans = spansByBound[bound].filter_sims(sim_i, each_n)
      var median_depth = medByBound[bound]

      if spans.len > 5_000:
        when defined(debug):
          stderr.write_line &"High depth for bound {opts.targets[bound.tid].name}:{bound.left}-{bound.right} got {spans.len} pairs. Skipping."
        continue
      if median_depth == -1:
        continue

      var gt = genotype(bound, str_reads, spans, opts, float(median_depth))

      var canon_repeat = bound.repeat.canonical_repeat
      if not genotypes_by_repeat.hasKey(canon_repeat):
        genotypes_by_repeat[canon_repeat] = @[]
      genotypes_by_repeat[canon_repeat].add(gt)

      #var estimate = spans.estimate_size(frag_dist)
      bounds_fh.write_line bound.tostring(opts.targets) & "\t" & $median_depth
      when defined(debug):
        for s in spans:
          span_fh.write_line s.tostring(bound, opts.targets[bound.tid].name)
        var locusid = bound.id(opts.targets)
        for r in str_reads:
          reads_fh.write_line r.tostring(opts.targets) & "\t" & $locusid

    # Loop through again and refine genotypes for loci that are the only
    # large expansion with that repeat unit
    for repeat, genotypes in genotypes_by_repeat:
      var gt_expanded: seq[Call]
      for gt in genotypes:
        if gt.is_large:
          gt_expanded.add(gt)
          if gt_expanded.len > 1:
            break
      if gt_expanded.len == 1:
        gt_expanded[0].update_genotype(unplaced_counts[repeat])
      for gt in genotypes:
        gt_fh.write_line gt.tostring()

    for repeat, count in unplaced_counts:
        unplaced_fh.write_line &"{repeat}\t{count}"

  gt_fh.close
  bounds_fh.close
  unplaced_fh.close
  when defined(debug):
    span_fh.close
    reads_fh.close
    stderr.write_line &"wrote str-like reads to {args.output_prefix}-reads.txt"
    stderr.write_line &"wrote spanning reads and spanning pairs to {args.output_prefix}-spanning.txt"
  if args.verbose:
    stderr.write_line &"Supporting evidence used to make the genotype calls:"
    stderr.write_line &"wrote putative str bounds to {args.output_prefix}-bounds.txt"
    stderr.write_line &"wrote counts of unplaced reads with STR content to {args.output_prefix}-unplaced.txt"
    stderr.write_line &"Main results file:"
    stderr.write_line &"wrote genotypes to {args.output_prefix}-genotype.txt"

when isMainModule:
  when not defined(danger):
   stderr.write_line "warning !!! not compiled in fast mode. Compile with -d:danger to increase speed"

  # Parse args/options
  call_main()

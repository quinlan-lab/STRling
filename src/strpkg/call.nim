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

proc `==`(a, b:Target): bool =
  return a.tid == b.tid and a.length == b.length and a.name == b.name

proc same(a:seq[Target], b:seq[Target]): bool =
  if a.len != b.len:
    return false
  for i, aa in a:
    if aa != b[i]:
      return false
  return true

proc call_main*() =
  var p = newParser("strling call"):
    option("-f", "--fasta", help="path to fasta file")
    option("-m", "--min-support", help="minimum number of supporting reads for a locus to be reported", default="6")
    option("-c", "--min-clip", help="minimum number of supporting clipped reads for each side of a locus", default="0")
    option("-t", "--min-clip-total", help="minimum total number of supporting clipped reads for a locus", default="0")
    option("-q", "--min-mapq", help="minimum mapping quality (does not apply to STR reads)", default="40")
    option("-l", "--loci", help="Annoated bed file specifying additional STR loci to genotype. Format is: chr start stop repeatunit [name]")
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
  # TODO: use other extracted stuffs here.
  for k, trs in treads_by_tid_rep.mpairs:
    trs.sort(proc(a, b:tread): int =
      cmp(a.position, b.position)
    )

  ### discovery
  var
    gt_fh:File
    reads_fh:File
    bounds_fh:File
    span_fh:File
    unplaced_fh:File
  if not open(gt_fh, args.output_prefix & "-genotype.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(reads_fh, args.output_prefix & "-reads.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(bounds_fh, args.output_prefix & "-bounds.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(span_fh, args.output_prefix & "-spanning.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(unplaced_fh, args.output_prefix & "-unplaced.txt", mode=fmWrite):
    quit "couldn't open output file"

  # Write headers
  bounds_fh.write_line(bounds_header & "\tdepth")
  reads_fh.write_line &"#chrom\tpos\tstr\tsoft_clip\tstr_count\tqname\tcluster_id"
  gt_fh.write_line("#chrom\tleft\tright\trepeatunit\tallele1_est\tallele2_est\toverlapping_reads\tspanning_reads\tspanning_pairs\tleft_clips\tright_clips\tunplaced_pairs\tdepth\tsum_str_counts")

  #TODO: read both a bounds file (from merge) and a bed file

  var loci: seq[Bounds]
  if args.loci != "":
    # Parse bed file of regions. These will also be genotyped
    loci = parse_bed(args.loci, opts.targets, opts.window.uint32)
    stderr.write_line &"Read {len(loci)} loci from {args.loci}"

  var unplaced_counts = initCountTable[string]()
  var genotypes_by_repeat = initTable[string, seq[Call]]()

  # Assign STR reads to provided loci and genotype them
  for locus in loci:
    var str_reads = assign_reads_locus(locus, treads_by_tid_rep)
    # Report spanning reads/pairs for any remaining loci that were not matched with a bound
    if locus.right - locus.left > 1000'u32:
      stderr.write_line "large bounds:" & $locus & " skipping"
      continue
    var (spans, median_depth) = ibam.spanners(locus, opts.window, frag_dist, opts.min_mapq)
    if spans.len > 5_000:
      when defined(debug):
        stderr.write_line &"High depth for bound {opts.targets[locus.tid].name}:{locus.left}-{locus.right} got {spans.len} pairs. Skipping."
      continue
    if median_depth == -1:
      continue

    var gt = genotype(locus, str_reads, spans, opts.targets, float(median_depth))

    var canon_repeat = locus.repeat.canonical_repeat
    if not genotypes_by_repeat.hasKey(canon_repeat):
      genotypes_by_repeat[canon_repeat] = @[]
    genotypes_by_repeat[canon_repeat].add(gt)

    #var estimate = spans.estimate_size(frag_dist)
    bounds_fh.write_line locus.tostring(opts.targets) & "\t" & $median_depth
    for s in spans:
      span_fh.write_line s.tostring(locus, opts.targets[locus.tid].name)


  # Cluster remaining reads and genotype
  var ci = 0

  for key, treads in treads_by_tid_rep.mpairs:
    ## treads are for the single tid, repeat unit defined by key
    for c in treads.cluster(max_dist=opts.window.uint32, min_supporting_reads=opts.min_support):
      if c.reads[0].tid == -1:
        unplaced_counts[c.reads[0].repeat.tostring] = c.reads.len
        continue

      var b: Bounds
      var good_cluster: bool
      (b, good_cluster) = check_cluster(c, min_clip, min_clip_total)
      if good_cluster == false:
        continue

      var (spans, median_depth) = ibam.spanners(b, opts.window, frag_dist, opts.min_mapq)
      if spans.len > 5_000:
        when defined(debug):
          stderr.write_line &"High depth for bound {opts.targets[b.tid].name}:{b.left}-{b.right} got {spans.len} pairs. Skipping."
        continue
      if median_depth == -1:
        continue

      var gt = genotype(b, c.reads, spans, opts.targets, float(median_depth))

      var canon_repeat = b.repeat.canonical_repeat
      if not genotypes_by_repeat.hasKey(canon_repeat):
        genotypes_by_repeat[canon_repeat] = @[]
      genotypes_by_repeat[canon_repeat].add(gt)

      #var estimate = spans.estimate_size(frag_dist)
      bounds_fh.write_line b.tostring(opts.targets) & "\t" & $median_depth

      for s in spans:
        span_fh.write_line s.tostring(b, opts.targets[b.tid].name)
      for s in c.reads:
        reads_fh.write_line s.tostring(opts.targets) & "\t" & $ci
      ci += 1

  # Loop through again and refine genotypes
  for repeat, genotypes in genotypes_by_repeat:
    if len(genotypes) == 1:
      var gt = genotypes[0]
      gt.update_genotype(unplaced_counts[repeat])
      gt_fh.write_line gt.tostring()
    else:
      for gt in genotypes:
        gt_fh.write_line gt.tostring()

  for repeat, count in unplaced_counts:
      unplaced_fh.write_line &"{repeat}\t{count}"

  gt_fh.close
  reads_fh.close
  bounds_fh.close
  span_fh.close
  unplaced_fh.close
  if args.verbose:
    stderr.write_line &"Supporting evidence used to make the genotype calls:"
    stderr.write_line &"wrote putative str bounds to {args.output_prefix}-bounds.txt"
    stderr.write_line &"wrote str-like reads to {args.output_prefix}-reads.txt"
    stderr.write_line &"wrote spanning reads and spanning pairs to {args.output_prefix}-spanning.txt"
    stderr.write_line &"wrote counts of unplaced reads with STR content to {args.output_prefix}-unplaced.txt"
    stderr.write_line &"Main results file:"
    stderr.write_line &"wrote genotypes to {args.output_prefix}-genotype.txt"

when isMainModule:
  when not defined(danger):
   stderr.write_line "warning !!! not compiled in fast mode. Compile with -d:danger to increase speed"

  # Parse args/options
  call_main()

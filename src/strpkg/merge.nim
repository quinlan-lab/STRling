import msgpack4nim
import hts/bam
import hts/fai
import strformat
import argparse
import algorithm

import ./cluster
import ./utils
import ./genome_strs
import ./extract
import ./callclusters
import ./unpack

export tread
export Soft

proc has_per_sample_reads(c:Cluster, supporting_reads:int, tread2sample:TableRef[tread, uint32]): bool =
  ## check that, within the cluster there are at least `supporting reads` from
  ## 1 sample.
  var sample_counts = initCountTable[uint32](8)
  for r in c.reads:
    sample_counts.inc(tread2sample[r])

  return sample_counts.largest().val >= supporting_reads


proc merge_main*() =
  var p = newParser("strling merge"):
    option("-f", "--fasta", help="path to fasta file (required if using CRAM input)")
    option("-w", "--window", help="Number of bp within which to search for reads supporting the other side of a bound. Estimated from the insert size distribution by default.", default = "-1")
    option("-m", "--min-support", help="minimum number of supporting reads required in at least one individual for a locus to be reported", default="6")
    option("-c", "--min-clip", help="minimum number of supporting clipped reads for each side of a locus", default="0")
    option("-t", "--min-clip-total", help="minimum total number of supporting clipped reads for a locus", default="0")
    option("-q", "--min-mapq", help="minimum mapping quality (does not apply to STR reads)", default="40")
    option("-l", "--bed", help="Annoated bed file specifying additional STR loci to genotype. Format is: chr start stop repeatunit [name]")
    option("-o", "--output-prefix", help="prefix for output files. Suffix will be -bounds.txt", default="strling")
    flag("-v", "--verbose")
    arg("bin", nargs = -1, help="One or more bin files previously created by `strling extract`")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "merge":
    argv = argv[1..argv.high]

  if len(argv) == 0: argv = @["-h"]

  var args = p.parse(argv)

  if args.help:
    quit 0

  var window = parseInt(args.window)
  var min_support = parseInt(args.min_support)
  var min_clip = uint16(parseInt(args.min_clip))
  var min_clip_total = uint16(parseInt(args.min_clip_total))
  var min_mapq = uint8(parseInt(args.min_mapq))
  var skip_reads = 100000

  if args.bed != "":
    if not fileExists(args.bed):
      quit "couldn't open bed file"

  var targets: seq[Target]
  var frag_dist: array[4096, uint32]
  var treads_by_tid_rep = newTable[tid_rep, seq[tread]](8192)
  # we use this to track (indirectly) which tread came from which sample.
  var tread2sample = newTable[tread, uint32]()

  for sample_i, binfile in args.bin:
    var fs = newFileStream(binfile, fmRead)
    var extracted = fs.unpack_file(drop_unplaced=true)
    fs.close()

    # TODO: Check all bin files came from the same version of STRling
    # TODO: Check all bin files used the same strling extract settings

    # Check all bin files are from the same reference genome (or at least have the same chroms)
    if targets.len == 0:
      targets = extracted.targets
    else:
      if not extracted.targets.same(targets):
        quit &"[strling] Error: inconsistent bam header for {binfile}. Were all samples run on the same reference genome?"

    # Aggregate insert sizes for all samples
    for i in 0..frag_dist.high:
      var before = frag_dist[i]
      frag_dist[i] += extracted.fragment_distribution[i]
      doAssert frag_dist[i] >= before, "overflow"

    # Unpack STR reads from all bin files and put them in a table by repeat unit and chromosome
    var sample_i = sample_i.uint32
    for r in extracted.reads.mitems:
      # set qname empty to save a bit of memory.
      r.qname = ""
      treads_by_tid_rep.mgetOrPut((r.tid, r.repeat), newSeqOfCap[tread](8192)).add(r)
      tread2sample[r] = sample_i
    #stderr.write_line &"[strling] read {extracted.reads.len} STR reads from file: {binfile}"
    for k, trs in treads_by_tid_rep.mpairs:
      trs.setLen(trs.len)

  var ntr = 0
  var n_unplaced = 0
  for k, trs in treads_by_tid_rep.mpairs:
    trs.sort(proc(a, b:tread): int =
      cmp(a.position, b.position)
    )
    ntr += trs.len
    for t in trs:
      if t.tid < 0: n_unplaced.inc
    trs.setLen(trs.len)
  stderr.write_line &"[strling] read {ntr} STR reads across all samples. unplaced: {n_unplaced} {100 * n_unplaced.float/ntr.float}%"

  if args.verbose:
    stderr.write_line "Calculated median fragment length accross all samples:", $frag_dist.median()
    stderr.write_line "10th, 90th percentile of fragment length:", $frag_dist.median(0.1), " ", $frag_dist.median(0.9)


  var opts = Options(median_fragment_length: frag_dist.median(0.98),
                      min_support: min_support, min_mapq: min_mapq, targets: targets)

  if window < 0:
    window = frag_dist.median(0.98)

  var loci: seq[Bounds]
  if args.bed != "":
    # Parse bed file of regions and report spanning reads
    loci = parse_bed(args.bed, targets, window.uint32)

  # Write bounds file with header
  var bounds_fh:File
  if not open(bounds_fh, args.output_prefix & "-bounds.txt", mode=fmWrite):
    quit "couldn't open output file"
  bounds_fh.write_line(bounds_header)

  # Assign STR reads to provided loci
  for locus in loci.mitems:
    var str_reads = assign_reads_locus(locus, treads_by_tid_rep)
    bounds_fh.write_line locus.tostring(opts.targets)

  # Cluster remaining reads
  var ci = 0
  for key, treads in treads_by_tid_rep.mpairs:
    shallow(treads)

    for c in treads.cluster(max_dist=window.uint32, min_supporting_reads=opts.min_support):
      if c.reads[0].tid == -1:
        continue
      if not c.has_per_sample_reads(opts.min_support, tread2sample): continue

      var b: Bounds
      var good_cluster: bool
      (b, good_cluster) = bounds(c, min_clip, min_clip_total)
      if good_cluster == false:
        continue

      bounds_fh.write_line b.tostring(targets)
      ci += 1

  bounds_fh.close
  if args.verbose:
    stderr.write_line &"Wrote merged str bounds to {args.output_prefix}-bounds.txt"

when isMainModule:
  when not defined(danger):
   stderr.write_line "warning !!! not compiled in fast mode. Compile with -d:danger to increase speed"

  # Parse args/options
  merge_main()

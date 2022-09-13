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

proc has_per_sample_reads(c:Cluster, supporting_reads:int): bool =
  ## check that, within the cluster there are at least `supporting reads` from
  ## 1 sample.
  var sample_counts = initCountTable[string](8)
  for r in c.reads:
    sample_counts.inc(r.qname)

  return sample_counts.largest().val >= supporting_reads

proc fill(targets:var seq[Target], fasta:string) =
  var fai:Fai
  if not fai.open(fasta):
    quit "could not open fasta:" & fasta

  for i in 0..<fai.len:
    let name = fai[i]
    targets.add(Target(tid:i.int, name:name, length:fai.chrom_len(name).uint32))


proc merge_main*() =
  var p = newParser("strling merge"):
    option("-f", "--fasta", help="path to fasta file (required if using CRAM input)")
    option("-w", "--window", help="Number of bp within which to search for reads supporting the other side of a bound. Estimated from the insert size distribution by default.", default = "-1")
    option("-m", "--min-support", help="minimum number of supporting reads required in at least one individual for a locus to be reported", default="5")
    option("-c", "--min-clip", help="minimum number of supporting clipped reads for each side of a locus", default="0")
    option("-t", "--min-clip-total", help="minimum total number of supporting clipped reads for a locus", default="0")
    option("-q", "--min-mapq", help="minimum mapping quality (does not apply to STR reads)", default="40")
    option("-l", "--bed", help="Annoated bed file specifying additional STR loci to genotype. Format is: chr start stop repeatunit [name]")
    option("-o", "--output-prefix", help="prefix for output files. Suffix will be -bounds.txt", default="strling")
    flag("-d", "--diff-refs", help="allow bin files generated on a mixture of reference genomes (by default differing references will produce an error). Reports chromosomes in the first bin or -f if provided")
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
  let allow_diff_chroms = args.diff_refs

  if args.bed != "":
    if not fileExists(args.bed):
      quit "couldn't open bed file"

  var targets: seq[Target]
  if args.fasta != "" and allow_diff_chroms:
    targets.fill(args.fasta)
  var frag_dist: array[4096, uint32]
  var treads_by_tid_rep = newTable[tid_rep, seq[tread]](8192)
  # we use this to track (indirectly) which tread came from which sample.

  for sample_i, binfile in args.bin:
    var fs = newFileStream(binfile, fmRead)
    if fs == nil:
      raise newException(IOError, &"[strling] unable to open {binFile} for reading. please check path")
    if args.verbose:
      stderr.write_line &"[strling] reading bin file: ", binfile
    var extracted = fs.unpack_file(drop_unplaced=true, verbose=args.verbose)
    fs.close()

    # Check all bin files are from the same reference genome (or at least have the same chroms)
    if targets.len == 0:
      targets = extracted.targets
    else:
      if not extracted.targets.same(targets) and not allow_diff_chroms:
        quit &"[strling] Error: inconsistent bam header for {binfile}. Were all samples run on the same reference genome?"

    # Aggregate insert sizes for all samples
    for i in 0..frag_dist.high:
      var before = frag_dist[i]
      frag_dist[i] += extracted.fragment_distribution[i]
      doAssert frag_dist[i] >= before, "overflow"

    # Unpack STR reads from all bin files and put them in a table by repeat unit and chromosome
    let sample_i = sample_i.uint32
    let sample_i_str = $sample_i

    for r in extracted.reads.mitems:
      # HACK: set qname to the sample_i so we can track reads per sample
      # this saves memory over having a separate tread->sample lookup.
      r.qname = sample_i_str
      treads_by_tid_rep.mgetOrPut((r.tid, r.repeat), newSeqOfCap[tread](128)).add(r)
      # NOTE: this doubles the memory because we have one entry for each read!
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
  if args.verbose:
    stderr.write_line &"[strling] read {ntr} STR reads across all samples."

  if args.verbose:
    stderr.write_line "[strling] Calculated median fragment length accross all samples:", $frag_dist.median()
    stderr.write_line "[strling] 10th, 90th percentile of fragment length:", $frag_dist.median(0.1), " ", $frag_dist.median(0.9)


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
      if not c.has_per_sample_reads(opts.min_support): continue

      var b: Bounds
      var good_cluster: bool
      let max_clip_dist: uint16 = uint16(0.5 * frag_dist.median(0.5).float)
      (b, good_cluster) = bounds(c, min_clip, min_clip_total, max_clip_dist)
      if good_cluster == false:
        continue

      bounds_fh.write_line b.tostring(targets)
      ci += 1

  bounds_fh.close
  if args.verbose:
    stderr.write_line &"[strling] Wrote merged str bounds to {args.output_prefix}-bounds.txt"

when isMainModule:
  when not defined(danger):
   stderr.write_line "[strling] WARNING!!! not compiled in fast mode. Compile with -d:danger to increase speed"

  # Parse args/options
  merge_main()

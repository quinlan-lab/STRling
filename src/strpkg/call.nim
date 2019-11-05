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
import ./cluster
import ./collect
import ./utils
import ./genotyper
import ./genome_strs
export tread
export Soft
import strformat
import math
import argparse
import ./extract

proc call_main*() =
  var p = newParser("strling call"):
    option("-f", "--fasta", help="path to fasta file")
    option("-m", "--min-support", help="minimum number of supporting reads for a locus to be reported", default="6")
    option("-c", "--min-clip", help="minimum number of supporting clipped reads for each side of a locus", default="2")
    option("-t", "--min-clip-total", help="minimum total number of supporting clipped reads for a locus", default="5")
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

  if not open(ibam_dist, args.bam, fai=args.fasta, threads=2):
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
  if not open(ibam, args.bam, fai=args.fasta, threads=2, index=true):
    quit "couldn't open bam"
  cram_opts = 8191 - SAM_RGAUX.int - SAM_QUAL.int
  discard ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, cram_opts)

  var cache = Cache(tbl:newTable[string, tread](8192), cache: newSeqOfCap[tread](65556))
  var opts = Options(median_fragment_length: frag_median,
                      min_support: min_support, min_mapq: min_mapq)
  var nreads = 0

  var fs = newFileStream(args.bin, fmRead)
  while not fs.atEnd:
    var t:tread
    fs.unpack(t)
    cache.cache.add(t)
  stderr.write_line &"[strling] read {cache.cache.len} treads from bin file"

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

  var window = frag_dist.median(0.98)

  reads_fh.write_line &"#chrom\tpos\tstr\tsoft_clip\tstr_count\tqname\tcluster_id" # print header
  gt_fh.write_line("#chrom\tleft\tright\trepeatunit\tallele1_est\tallele2_est\ttotal_reads\tspanning_reads\tspanning_pairs\tleft_clips\tright_clips\tunplaced_pairs\tdepth\tsum_str_counts")



  var targets = ibam.hdr.targets

  var loci: seq[Bounds]
  if args.loci != "":
    # Parse bed file of regions and report spanning reads
    loci = parse_loci(args.loci, targets)

  var unplaced_counts = initCountTable[string]()
  var genotypes_by_repeat = initTable[string, seq[Call]]()

  var ci = 0
  for c in cache.cache.cluster(max_dist=window.uint32, min_supporting_reads=opts.min_support):
    if c.reads[0].tid == -1:
      unplaced_fh.write_line &"{c.reads[0].repeat.tostring}\t{c.reads.len}"
      unplaced_counts[c.reads[0].repeat.tostring] = c.reads.len
      continue
    if c.reads.len >= uint16.high.int:
      stderr.write_line "More than " & &"{uint16.high.int}" & " reads in cluster with first read:" & $c.reads[0] & " skipping"
      continue
    var b = c.bounds

    # Check if bounds overlaps with one of the input loci, if so overwrite b attributes
    for i, locus in loci:
      if b.overlaps(locus):
        b.name = locus.name
        b.left = locus.left
        b.right = locus.right
        # Remove locus from loci (therefore will use first matching bound and locus)
        loci.del(i)
        break

    if b.right - b.left > 1000'u32:
      stderr.write_line "large bounds:" & $b & " skipping"
      continue
    # require left and right support
    if b.n_left < min_clip: continue
    if b.n_right < min_clip: continue
    if (b.n_right + b.n_left) < min_clip_total: continue

    #[
    #debugging for reference size
    var ref_size = b.check_reference_size(fai, targets[b.tid].name)
    echo ref_size, b

    var gr: seq[region]

    try:
      if genome_str[targets[b.tid].name].find(b.left.int - 10, b.right.int + 10, gr):
      else:
    except KeyError:
      discard
    ]#

    var (spans, median_depth) = ibam.spanners(b, window, frag_dist, opts.min_mapq)
    var gt = genotype(b, c.reads, spans, targets, float(median_depth))

    var canon_repeat = b.repeat.canonical_repeat
    if not genotypes_by_repeat.hasKey(canon_repeat):
      genotypes_by_repeat[canon_repeat] = @[]
    genotypes_by_repeat[canon_repeat].add(gt)

    var estimate = spans.estimate_size(frag_dist)
    bounds_fh.write_line b.tostring(targets) & "\t" & $median_depth
    for s in spans:
      span_fh.write_line s.tostring(b, targets[b.tid].name)
    for s in c.reads:
      reads_fh.write_line s.tostring(targets) & "\t" & $ci
    ci += 1

  # Report spanning reads/pairs for any remaining loci that were not matched with a bound
  for locus in loci:
    if locus.right - locus.left > 1000'u32:
      stderr.write_line "large bounds:" & $locus & " skipping"
      continue
    var (spans, median_depth) = ibam.spanners(locus, window, frag_dist, opts.min_mapq)
    var estimate = spans.estimate_size(frag_dist)
    bounds_fh.write_line locus.tostring(targets) & "\t" & $median_depth
    for s in spans:
      span_fh.write_line s.tostring(locus, targets[locus.tid].name)
 
  ### end discovery

  # Loop through again and refine genotypes
  for repeat, genotypes in genotypes_by_repeat:
    if len(genotypes) == 1:
      var gt = genotypes[0]
      gt.update_genotype(unplaced_counts[repeat])
      gt_fh.write_line gt.tostring()
    else:
      for gt in genotypes:
        gt_fh.write_line gt.tostring()

  gt_fh.close
  reads_fh.close
  bounds_fh.close
  span_fh.close
  unplaced_fh.close
  if args.verbose:
    stderr.write_line cache.tbl.len, " left in table"
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

  call_main()


  # Parse args/options

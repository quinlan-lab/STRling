import msgpack4nim
import hts/bam
import hts/fai
import ./cluster
import ./utils
import ./genome_strs
import strformat
import argparse
import ./extract
export tread
export Soft

proc merge_main*() =
  var p = newParser("strling merge"):
    option("-f", "--fasta", help="path to fasta file (required if using CRAM input)")
    option("-w", "--window", help="Number of bp within which to search for reads supporting the other side of a bound. Estimated from the insert size distribution by default.", default = "-1")
    option("-m", "--min-support", help="minimum number of supporting reads for a locus to be reported", default="6")
    option("-c", "--min-clip", help="minimum number of supporting clipped reads for each side of a locus", default="0")
    option("-t", "--min-clip-total", help="minimum total number of supporting clipped reads for a locus", default="0")
    option("-q", "--min-mapq", help="minimum mapping quality (does not apply to STR reads)", default="40")
    option("-l", "--loci", help="Annoated bed file specifying additional STR loci to genotype. Format is: chr start stop repeatunit [name]")
    option("-o", "--output-prefix", help="prefix for output files", default="strling")
    flag("-v", "--verbose")
    arg("bam", help="path to a single representitive bam or cram file (for extracting fragment size dist)")
    arg("bin", nargs = -1, help="One or more bin files previously created by `strling extract`")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "merge":
    argv = argv[1..argv.high]

  if len(argv) == 0: argv = @["-h"]

  var args = p.parse(argv)

  if args.help:
    quit 0

  var ibam_dist:Bam
  var window = parseInt(args.window)
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
  var frag_median = frag_dist.median
  var targets = ibam_dist.hdr.targets
  if args.verbose:
    stderr.write_line "Calculated median fragment length:", frag_median
    stderr.write_line "10th, 90th percentile of fragment length:", $frag_dist.median(0.1), " ", $frag_dist.median(0.9)
  ibam_dist.close()
  ibam_dist = nil

  var cache = Cache(tbl:newTable[string, tread](8192), cache: newSeqOfCap[tread](65556))
  var opts = Options(median_fragment_length: frag_median,
                      min_support: min_support, min_mapq: min_mapq)

  for binfile in args.bin:
    var fs = newFileStream(binfile, fmRead)
    var nreads = 0
    while not fs.atEnd:
      var t:tread
      fs.unpack(t)
      cache.cache.add(t)
      nreads += 1
    stderr.write_line &"[strling] read {nreads} STR reads from file: {binfile}"

  var bounds_fh:File
  if not open(bounds_fh, args.output_prefix & "-bounds.txt", mode=fmWrite):
    quit "couldn't open output file"

  #XXX this should already have been estimated for all samples? Record in bin and take average?
  if window < 0:
    window = frag_dist.median(0.98)

  var loci: seq[Bounds]
  if args.loci != "":
    # Parse bed file of regions and report spanning reads
    loci = parse_loci(args.loci, targets)

  # Write header
  bounds_fh.write_line "#chrom\tleft\tright\trepeat\tname\tcenter_mass\tn_left\tn_right\tn_total"

  var ci = 0
  for c in cache.cache.cluster(max_dist=window.uint32, min_supporting_reads=opts.min_support):
    if c.reads[0].tid == -1:
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
        b.force_report = true
        # Remove locus from loci (therefore will use first matching bound and locus)
        loci.del(i)
        break

    if b.right - b.left > 1000'u32:
      stderr.write_line "large bounds:" & $b & " skipping"
      continue
    # require left and right support
    if not b.force_report:
      if b.n_left < min_clip: continue
      if b.n_right < min_clip: continue
      if (b.n_right + b.n_left) < min_clip_total: continue

    bounds_fh.write_line b.tostring(targets)
    ci += 1

  bounds_fh.close
  if args.verbose:
    stderr.write_line cache.tbl.len, " left in table"
    stderr.write_line &"Wrote merged str bounds to {args.output_prefix}-bounds.txt"

when isMainModule:
  when not defined(danger):
   stderr.write_line "warning !!! not compiled in fast mode. Compile with -d:danger to increase speed"

  # Parse args/options
  merge_main()

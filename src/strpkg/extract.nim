import kmer
import algorithm
import ./genome_strs
import msgpack4nim
import strutils
import times
import random
import tables
import hts/bam
import ./cluster
import ./collect
import ./utils
export tread
export Soft
import strformat
import math
import argparse

proc get_repeat*(aln:Record, genome_str:TableRef[string, Lapper[region]], counts:var Seqs[uint8], repeat_count: var int, align_length: var int, opts:Options): array[6, char] =
  repeat_count = 0
  when defined(skipFullMatch):
    # compile with -d:skipFullMatch to make it much faster
    if aln.cigar.len == 1 and aln.cigar[0].op == CigarOp.match:
      align_length = aln.cigar[0].len
      return

  # don't check for strs when we know there's not repeat in the reference *and*
  # we have an exact match to the reference.
  if aln.cigar.len == 1 and aln.cigar[0].op == CigarOp.match and aln.chrom in genome_str:
    var empty: seq[region]
    if aln.chrom notin genome_str or not genome_str[aln.chrom].find(aln.start, aln.stop, empty):
      align_length = aln.cigar[0].len
      return

  var read = ""
  aln.sequence(read)
  align_length = len(read)

  if aln.cigar.len > 0:
    # we only test the aligned part of the read for repeats.
    # if it is soft-clipped, those are checked separately anyway.
    if aln.cigar[0].op == CigarOp.soft_clip:
      read = read[aln.cigar[0].len..<read.len]
      align_length -= aln.cigar[0].len

    var L = aln.cigar.len
    if aln.cigar[L-1].op == CigarOp.soft_clip:
      read = read[0..<read.len-aln.cigar[L-1].len]
      align_length -= aln.cigar[L-1].len


  result = read.get_repeat(counts, repeat_count, opts, false)


proc tostring*(t:tread, targets: seq[Target]): string =
  var chrom = if t.tid == -1: "unknown" else: targets[t.tid].name
  var rep: string
  for v in t.repeat:
    if v == 0.char: continue
    rep.add(v)
  result = &"""{chrom}	{t.position}	{rep}	{t.split}	{t.repeat_count}	{t.qname}"""

proc repeat_length(t:tread): uint8 {.inline.} =
  for v in t.repeat:
    if v == 0.char: return
    result.inc

template p_repeat(t:tread): float =
  # proportion repeat
  float(t.repeat_count * t.repeat_length) / t.align_length.float

template after_mate(aln:Record): bool {.dirty.} =
  (aln.tid > aln.mate_tid or (aln.tid == aln.mate_tid and ((aln.start > aln.mate_pos) or (aln.start == aln.mate_pos and cache.tbl.hasKey(aln.qname)))))

proc to_tread(aln:Record, genome_str:TableRef[string, Lapper[region]], counts: var Seqs[uint8], opts:Options): tread {.inline.} =
  var repeat_count: int
  var align_length: int
  var rep = aln.get_repeat(genome_str, counts, repeat_count, align_length, opts)
  doAssert align_length > 0, aln.tostring
  doAssert repeat_count < 256, aln.tostring

  result = tread(tid:aln.tid.int32,
                 position: max(0, aln.start).uint32,
                 repeat: rep,
                 flag: aln.flag,
                 repeat_count: repeat_count.uint8,
                 align_length: align_length.uint8,
                 split: Soft.none,
                 mapping_quality: aln.mapping_quality,
                 qname: aln.qname)

type Cache* = object
  tbl*: TableRef[string, tread]
  cache*: seq[tread]

proc add_soft*(cache:var Cache, aln:Record, counts: var Seqs[uint8], opts:Options, read_repeat:array[6, char]) =
  # if read_repeat is not empty, we are more permissive on the length of the
  # soft-clip as long as it has the same repeat unit.

  if aln.mapping_quality < opts.min_mapq: return
  if aln.cigar.len == 0 or (aln.cigar[0].op != CigarOp.soft_clip  and aln.cigar[aln.cigar.len - 1].op != CigarOp.soft_clip): return

  var soft_seq: string
  var repeat_count: int
  for cig_index in [0, aln.cigar.len - 1]:
    var c = aln.cigar[cig_index]
    if c.op != CigarOp.soft_clip: continue

    if read_repeat[0] == 0.char and c.len <= 16: continue

    aln.sequence(soft_seq)
    if cig_index == 0:
      soft_seq = soft_seq[0..<c.len]
    else:
      soft_seq = soft_seq[soft_seq.len - c.len..<soft_seq.len]

    var repeat = soft_seq.get_repeat(counts, repeat_count, opts)
    if repeat_count == 0: continue

    # If read is soft-clipped on the left take the read position as the start of read
    # If soft-clipped on the right take the read position to the the end of the alignment
    var position = if cig_index == 0: max(0, aln.start).uint32 else: max(0, aln.stop).uint32
    cache.cache.add(tread(tid:aln.tid.int32,
                  position: position,
                  flag: aln.flag,
                  repeat: repeat,
                  repeat_count: repeat_count.uint8,
                  align_length: soft_seq.len.uint8,
                  split: if cig_index == 0: Soft.left else: Soft.right,
                  mapping_quality: aln.mapping_quality,
                  qname: aln.qname
                  ))

proc should_reverse(f:Flag): bool {.inline.} =
  ## this is only  called after we've ensured hi-quality of mate and lo-quality
  ## of self.
  result = not f.mate_reverse
  if f.reverse:
    result = not result

proc adjust_by(A:var tread, B:tread, opts:Options): bool =
  if A.repeat_count == 0'u8: return false
  # potentially adjust A position by B

  result = true

  # when B has hi mapping quality, we adjust A if:
  # A is very repetitive and B is not very repetitive
  # A is mapped poorly, B is mapped well and it's not a proper pair
  # TODO: use opts.$param for 0.7 and 0.2
  if B.mapping_quality > opts.min_mapq and ((A.p_repeat > 0.7 and B.p_repeat < 0.2) or (not A.flag.proper_pair and A.mapping_quality < opts.min_mapq)):
    # B is right of A, so the position subtracts th fragment length
    # Note fragment size is the external distance
    if B.flag.reverse:
      #   A                  B
      #   =======>             <===========
      #   000000000000000000000000000000000 fragment length
      A.position = B.position - opts.median_fragment_length.uint32 + B.align_length + uint32(A.align_length.float / 2'f + 0.5)
    else:
      #   B                 A
      #   =======>             <===========
      #   000000000000000000000000000000000 fragment length
      A.position = B.position + opts.median_fragment_length.uint32 - uint32(A.align_length.float / 2'f + 0.5)

    A.tid = B.tid
    A.mapping_quality = max(A.mapping_quality, B.mapping_quality)
    if A.flag.should_reverse:
      A.repeat.min_rev_complement

  # A is STR, A is mapped well
  elif A.mapping_quality >= opts.min_mapq or A.flag.proper_pair:
    A.position += uint32(A.align_length.float / 2'f + 0.5) # Record position as middle of A
    A.mapping_quality = max(A.mapping_quality, B.mapping_quality)

# Return true if both reads in pair are STR, or one is STR, one low mapping qual
proc unplaced_pair*(A:var tread, B:tread, opts:Options): bool =
  if A.p_repeat > opts.proportion_repeat and B.p_repeat > opts.proportion_repeat:
    return true
  if A.p_repeat > opts.proportion_repeat and B.mapping_quality < opts.min_mapq:
    return true
  if B.p_repeat > opts.proportion_repeat and A.mapping_quality < opts.min_mapq:
    return true

  return false

proc add(cache:var Cache, aln:Record, genome_str:TableRef[string, Lapper[region]], counts: var Seqs[uint8], opts:var Options) =
  doAssert not (aln.flag.secondary or aln.flag.supplementary)

  # Check if you have both reads
  if aln.after_mate:
    var mate:tread
    if not cache.tbl.take(aln.qname, mate): return

    # get mates together and see what, if anything should be added to the final
    # cache.

    # set the aln qual so that it transfers to soft and self
    var self = aln.to_tread(genome_str, counts, opts)
    # drop the proportion_repeat for the soft-clipped portion as it
    # often has imperfect repeats.
    var b = opts.proportion_repeat
    opts.proportion_repeat = min(b, 0.6)
    cache.add_soft(aln, counts, opts, self.repeat)
    # then set it back
    opts.proportion_repeat = b
    if mate.repeat_count == 0'u8 and self.repeat_count == 0: return

    # If both reads in pair are STR, or one is STR, one low mapping qual, set position to unknown
    if unplaced_pair(self, mate, opts):
      # we only keep unplaced pairs if both reads have some repeat unit
      # we could require them to have the same repeat unit...
      if self.repeat[0] == '\0' or mate.repeat[0] == '\0':
        return

      # Since the true location of these reads is unknown, set the 
      # repeat to the first alphabetically of seq and reverse complement
      self.repeat = self.repeat.canonical_repeat
      self.position = uint32(0)
      self.tid = -1
      mate.repeat = mate.repeat.canonical_repeat
      mate.position = uint32(0)
      mate.tid = -1
      cache.cache.add(self)
      cache.cache.add(mate)
      return

    if mate.adjust_by(self, opts):
      cache.cache.add(mate)
    if self.adjust_by(mate, opts):
      cache.cache.add(self)

  else:
    var tr = aln.to_tread(genome_str, counts, opts)
    var b = opts.proportion_repeat
    opts.proportion_repeat -= 0.07
    cache.add_soft(aln, counts, opts, tr.repeat)
    opts.proportion_repeat = b
    if cache.tbl.hasKeyOrPut(aln.qname, tr):
      stderr.write_line "[strling] warning. bad read (this happens with bwa-kit alignments):" & aln.qname & " already in table as:" & $cache.tbl[aln.qname]
      var mate:tread
      discard cache.tbl.take(aln.qname, mate)

proc extract_main*() =
  # Parse args/options
  var p = newParser("str extract"):
    option("-f", "--fasta", help="path to fasta file (required for CRAM)")
    option("-g", "--genome-repeats", help="optional path to genome repeats file. if it does not exist, it will be created")
    option("-p", "--proportion-repeat", help="proportion of read that is repetitive to be considered as STR", default="0.8")
    option("-q", "--min-mapq", help="minimum mapping quality (does not apply to STR reads)", default="40")
    flag("-v", "--verbose")
    arg("bam", help="path to bam file")
    arg("bin", help="path bin to output bin file to be created")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "extract":
    argv = argv[1..argv.high]
  if len(argv) == 0: argv = @["--help"]

  var args = p.parse(argv)
  if args.help:
    quit 0

  var ibam:Bam
  var proportion_repeat = parseFloat(args.proportion_repeat)
  var min_mapq = uint8(parseInt(args.min_mapq))
  var skip_reads = 100000

  if not open(ibam, args.bam, fai=args.fasta, threads=2):
    quit "couldn't open bam"

  var cram_opts = 8191 - SAM_RNAME.int - SAM_RGAUX.int - SAM_QUAL.int - SAM_SEQ.int
  discard ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, cram_opts)

  var frag_dist = ibam.fragment_length_distribution(skip_reads=skip_reads)
  #echo frag_dist[0..<1000]
  var frag_median = frag_dist.median
  if args.verbose:
    stderr.write_line "Calculated median fragment length:", frag_median
    stderr.write_line "10th, 90th percentile of fragment length:", $frag_dist.median(0.1), " ", $frag_dist.median(0.9)

  ibam.close()
  if not open(ibam, args.bam, fai=args.fasta, threads=2, index=true):
    quit "couldn't open bam"
  cram_opts = 8191 - SAM_RGAUX.int - SAM_QUAL.int
  discard ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, cram_opts)

  var decodeds = newSeq[string](7)
  for i, s in decodeds.mpairs:
    decodeds[i] = newString(i)
  shallow(decodeds)

  var cache = Cache(tbl:newTable[string, tread](8192), cache: newSeqOfCap[tread](65556))
  var opts = Options(median_fragment_length: frag_median, proportion_repeat: proportion_repeat,
                      min_mapq: min_mapq)
  var nreads = 0
  var genome_str = genome_repeats(args.fasta, opts, args.genome_repeats)

  var t0 = cpuTime()
  var counts = init[uint8]()
  stderr.write_line "[strling] collecting str-like reads"
  var tid = -1
  for aln in ibam: #.query("14:92537254-92537477"):
    if aln.flag.secondary or aln.flag.supplementary: continue
    if aln.tid != tid and aln.tid >= 0:
      if ibam.hdr.targets[aln.tid].length > 2_000_000'u32:
        stderr.write_line "[strling] extracting chromosome:", $aln.chrom
      tid = aln.tid

    nreads.inc

    if nreads mod 10_000_000 == 0:
      var nrps = nreads.float64 / (cpuTime() - t0)
      if args.verbose:
        stderr.write_line $nreads, &" @ {aln.chrom}:{aln.start} {nrps:.1f} reads/sec tbl len: {cache.tbl.len} cache len: {cache.cache.len}"

    cache.add(aln, genome_str, counts, opts)

  stderr.write_line "[strling] extracting unampped reads"
  # get unmapped reads
  for aln in ibam.query("*"):
    if aln.flag.secondary or aln.flag.supplementary: continue
    nreads.inc
    cache.add(aln, genome_str, counts, opts)

  stderr.write_line "[strling] writing binary file:", args.bin
  var fs = newFileStream(args.bin, fmWrite)
  if fs == nil:
    quit "[strling] couldnt open binary output file"
  # TODO: write min_mapq, proportion repeat to start of bin file
  for c in cache.cache:
    fs.pack(c)
  stderr.write_line "[strling] finished extraction"
  fs.close

  ibam.close

when isMainModule:
  extract_main()

import kmer
import algorithm
import strutils
import times
import random
import tables
import hts/bam
import ./strpkg/cluster
import ./strpkg/collect
import ./strpkg/utils
export tread
export Soft
import strformat
import math
import argparse

# This is the bottleneck for run time at the moment
proc get_repeat(read: var string, counts: var Seqs[uint8], repeat_count: var int, opts:Options): array[6, char] =
  repeat_count = 0
  if read.count('N') > 20: return
  var s = newString(6)

  var best_score: int = -1
  for k in 2..6:
    let count = read.count(k, counts[k])
    var score = count * k
    if score <= best_score:
      if count < (read.len.float * 0.12 / k.float).int:
        break
      continue
    best_score = score
    if count > (read.len.float * opts.proportion_repeat / k.float).int:
      # now check the actual string because the kmer method can't track phase
      s = newString(k)
      counts[k].argmax.decode(s)
      if read.count(s) > (read.len.float * opts.proportion_repeat / k.float).int:
        copyMem(result.addr, s[0].addr, k)
        repeat_count = count
    elif count < (read.len.float * 0.12 / k.float).int:
      # e.g. for a 5 mer repeat, we should see some 2, 3, 4 mer reps and we can
      # bail if we do not. this is an optimization to avoid counting when we
      # cant possibly see a repeat.
      break

proc get_repeat*(aln:Record, counts:var Seqs[uint8], repeat_count: var int, align_length: var int, opts:Options): array[6, char] =
  # returns blank array if nothing passes.
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

  result = read.get_repeat(counts, repeat_count, opts)

proc tostring*(t:tread, targets: seq[Target]): string =
  var chrom = if t.tid == -1: "unknown" else: targets[t.tid].name
  var rep: string
  for v in t.repeat:
    if v == 0.char: continue
    rep.add(v)
  result = &"""{chrom}	{t.position}	{rep}	{t.split}	{t.repeat_count}"""
  when defined(debug) or defined(qname):
    result &= "\t" & t.qname

proc repeat_length(t:tread): uint8 {.inline.} =
  for v in t.repeat:
    if v == 0.char: return
    result.inc

template p_repeat(t:tread): float =
  # proportion repeat
  float(t.repeat_count * t.repeat_length) / t.align_length.float

template after_mate(aln:Record): bool {.dirty.} =
  (aln.tid > aln.mate_tid or (aln.tid == aln.mate_tid and ((aln.start > aln.mate_pos) or (aln.start == aln.mate_pos and cache.tbl.hasKey(aln.qname)))))

proc to_tread(aln:Record, counts: var Seqs[uint8], opts:Options): tread {.inline.} =
  var repeat_count: int
  var align_length: int
  var rep = aln.get_repeat(counts, repeat_count, align_length, opts)
  doAssert align_length > 0, aln.tostring
  doAssert repeat_count < 256, aln.tostring

  result = tread(tid:aln.tid.int32,
                 position: max(0, aln.start).uint32,
                 repeat: rep,
                 flag: aln.flag,
                 repeat_count: repeat_count.uint8,
                 align_length: align_length.uint8,
                 split: Soft.none,
                 mapping_quality: aln.mapping_quality)
  when defined(debug) or defined(qname):
    result.qname = aln.qname

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

    if read_repeat[0] == 0.char and c.len < 20: continue

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
                  mapping_quality: aln.mapping_quality
                  ))
    when defined(debug) or defined(qname):
      cache.cache[cache.cache.high].qname = aln.qname

proc should_reverse(f:Flag): bool {.inline.} =
  ## this is only  called after we've ensured hi-quality of mate and lo-quality
  ## of self.
  result = not f.mate_reverse
  if f.reverse:
    result = not result

proc min_rev_complement(repeat: var array[6, char]) {.inline.} =
  # find the minimal reverse complement of this repeat.
  # this turns the repeat in array to a string, doubles it so that ACT becomes
  # (the complement of) ACTACT and finds the minimal 3-mer of ACT, CTA, TAC
  # this could be optimized but probably won't be called often
  var s: string
  for c in repeat:
    if c == 0.char: break
    s.add(c)
  s = s.reverse_complement
  let l = s.len
  s.add(s)
  var mv = uint64(0) - 1'u64
  for m in s.slide_by(l):
    if m < mv:
      mv = m
  var ms = newString(l)
  mv.decode(ms)
  for i in 0..<l:
    repeat[i] = ms[i]

proc adjust_by(A:var tread, B:tread, opts:Options): bool =
  if A.repeat_count == 0'u8: return false
  # potentially adjust A position by B

  result = true

  # when B has hi mapping quality, we adjust A if:
  # A is very repetitive and B is not very repetitive
  # A is mapped poorly, B is mapped well and it's not a proper pair
  # TODO: use opts.$param for 0.9 and 0.2
  if B.mapping_quality > opts.min_mapq and ((A.p_repeat > 0.9 and B.p_repeat < 0.2) or (not A.flag.proper_pair and A.mapping_quality < opts.min_mapq)):
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
    if A.flag.should_reverse:
      A.repeat.min_rev_complement

  # A is STR, A is mapped well
  elif A.mapping_quality >= opts.min_mapq or A.flag.proper_pair:
    A.position += uint32(A.align_length.float / 2'f + 0.5) # Record position as middle of A

# Return true if both reads in pair are STR, or one is STR, one low mapping qual
proc unplaced_pair*(A:var tread, B:tread, opts:Options): bool =
  if A.p_repeat > opts.proportion_repeat and B.p_repeat > opts.proportion_repeat:
    return true
  if A.p_repeat > opts.proportion_repeat and B.mapping_quality < opts.min_mapq:
    return true
  if B.p_repeat > opts.proportion_repeat and A.mapping_quality < opts.min_mapq:
    return true

  return false

proc add(cache:var Cache, aln:Record, counts: var Seqs[uint8], opts:Options) =
  doAssert not (aln.flag.secondary or aln.flag.supplementary)

  # Check if you have both reads
  if aln.after_mate:
    var mate:tread
    if not cache.tbl.take(aln.qname, mate): return

    # get mates together and see what, if anything should be added to the final
    # cache.

    var self = aln.to_tread(counts, opts)
    cache.add_soft(aln, counts, opts, self.repeat)
    if mate.repeat_count == 0'u8 and self.repeat_count == 0: return

    # If both reads in pair are STR, or one is STR, one low mapping qual, set position to unknown
    if unplaced_pair(self, mate, opts):
      # NOTE: we don't know if we need to reverse complement these reads
      # so we will have to equate forward and reverse repeat units later.
      self.position = uint32(0)
      self.tid = -1
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
    var tr = aln.to_tread(counts, opts)
    cache.add_soft(aln, counts, opts, tr.repeat)
    doAssert not cache.tbl.hasKeyOrPut(aln.qname, tr), "error with read:" & aln.qname & " already in table as:" & $cache.tbl[aln.qname]


proc tostring*(a:array[6, char]): string =
  for c in a:
    if c == 0.char: return
    result.add(c)

when isMainModule:
  import math

  when not defined(danger):
   stderr.write_line "warning !!! not compiled in fast mode. Compile with -d:danger to increase speed"

  #testing()
  #if true:
  #  quit "asdf"
  #
  when defined(debug):
    block:
      var cnt = Seq[uint8](A: newSeq[uint8](4096))
      # 5mer
      for read in ["ATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCGATCCG",
      # 4mer
      "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT",
      "AGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACATGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGAC",
      "AGACAGATAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGAC",
      ]:
        var read = read
        stderr.write_line "read:", read
        #shuffle(read)
        for k in 2..6:
          cnt.clear
          var s = newString(k)
          var count = read.count(k, cnt)
          cnt.argmax.decode(s)
          stderr.write_line "k:", k, " ", count, " ", s, " break?:", count < (read.len.float * 0.15 / k.float).int, " imax:", cnt.imax, " cutoff:", int(read.len.float * 0.15 / k.float), " score:", k * count

  # Parse args/options
  var p = newParser("str"):
    option("-f", "--fasta", help="path to fasta file")
    option("-p", "--proportion-repeat", help="proportion of read that is repetitive to be considered as STR", default="0.8")
    option("-m", "--min-support", help="minimum number of supporting reads for a locus to be reported", default="5")
    option("-q", "--min-mapq", help="minimum mapping quality (does not apply to STR reads)", default="20")
    option("--skip", "Skip this many reads before calculating the insert size distribution", default="100000")
    option("-o", "--output-prefix", help="prefix for output files", default="strstrstr")
    flag("-v", "--verbose")
    arg("bam", help="path to bam file")

  var argv = commandLineParams()
  if len(argv) == 0: argv = @["-h"]
  var args = p.parse(argv)
  if args.help:
    quit 0

  var t0 = cpuTime()
  var ibam_dist:Bam
  var proportion_repeat = parseFloat(args.proportion_repeat)
  var min_support = parseInt(args.min_support)
  var min_mapq = uint8(parseInt(args.min_mapq))
  var skip_reads = parseInt(args.skip)

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


  var decodeds = newSeq[string](7)
  for i, s in decodeds.mpairs:
    decodeds[i] = newString(i)
  shallow(decodeds)

  var cache = Cache(tbl:newTable[string, tread](8192), cache: newSeqOfCap[tread](65556))
  var opts = Options(median_fragment_length: frag_median, proportion_repeat: proportion_repeat,
                      min_support: min_support, min_mapq: min_mapq)

  var nreads = 0
  var counts = init[uint8]()
  stderr.write_line "collecting str-like reads"
  for aln in ibam: #.query("14:92537254-92537477"):
    if aln.flag.secondary or aln.flag.supplementary: continue

    nreads.inc

    if nreads mod 10_000_000 == 0:
      var nrps = nreads.float64 / (cpuTime() - t0)
      if args.verbose:
        stderr.write_line $nreads, &" @ {aln.chrom}:{aln.start} {nrps:.1f} reads/sec tbl len: {cache.tbl.len} cache len: {cache.cache.len}"

    cache.add(aln, counts, opts)

  # get unmapped reads
  for aln in ibam.query("*"):
    if aln.flag.secondary or aln.flag.supplementary: continue
    nreads.inc
    cache.add(aln, counts, opts)

  stderr.write_line "[str] done reading bam, starting clustering"
  var
    reads_fh:File
    bounds_fh:File
    span_fh:File
    unplaced_fh:File
  if not open(reads_fh, args.output_prefix & "-reads.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(bounds_fh, args.output_prefix & "-bounds.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(span_fh, args.output_prefix & "-spanning.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(unplaced_fh, args.output_prefix & "-unplaced.txt", mode=fmWrite):
    quit "couldn't open output file"

  var window = frag_dist.median(0.98)

  reads_fh.write_line "chrom\tpos\tstr\tsoft_clip\tstr_count\tqname\tcluster_id" # print header
  var targets = ibam.hdr.targets
  var ci = 0
  for c in cache.cache.cluster(max_dist=window.uint32, min_supporting_reads=opts.min_support):
    if c.reads[0].tid == -1:
      unplaced_fh.write_line &"{c.reads[0].repeat.tostring}\t{c.reads.len}"
      continue
    var b = c.bounds
    if b.right - b.left > 1000'u32:
      stderr.write_line "large bounds:" & $b & " skipping"
      continue
    var spans = ibam.spanners(b, window, frag_dist, opts.min_mapq)
    var estimate = spans.estimate_size(frag_dist)
    bounds_fh.write_line b.tostring(targets) & "\testimate:" & $estimate
    for s in spans:
      span_fh.write_line s.tostring(b, targets[b.tid].name)
    for s in c.reads:
      reads_fh.write_line s.tostring(targets) & "\t" & $ci
    ci += 1

  reads_fh.close
  bounds_fh.close
  span_fh.close
  unplaced_fh.close
  if args.verbose:
    stderr.write_line cache.tbl.len, " left in table"
    stderr.write_line &"wrote bounds to {args.output_prefix}-bounds.txt"
    stderr.write_line &"wrote reads to {args.output_prefix}-reads.txt"
    stderr.write_line &"wrote spanners to {args.output_prefix}-spanning.txt"
    stderr.write_line &"wrote counts of unplaced fragments with STR content to {args.output_prefix}-unplaced.txt"

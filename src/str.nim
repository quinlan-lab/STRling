import kmer
import algorithm
import strutils
import times
import random
import tables
import hts/bam
import ./strpkg/cluster
import ./strpkg/utils
export tread
import strformat
import math
import argparse

type Seq[T] = object
  imax: int
  A: seq[T]

type Seqs[T] = array[7, Seq[T]]

type Options = object
  median_fragment_length: int
  proportion_repeat: float
  min_mapq: uint8


proc init[T](): Seqs[T] =
  result = [
     Seq[T](A: newSeq[T](0)),
     Seq[T](A: newSeq[T](0)),
     Seq[T](A: newSeq[T](16)),
     Seq[T](A: newSeq[T](64)),
     Seq[T](A: newSeq[T](256)),
     Seq[T](A: newSeq[T](1024)),
     Seq[T](A: newSeq[T](4096)),
  ]

proc fragment_length_distribution(bam:Bam, n_reads:int=2_000_000, skip_reads:int=100_000): array[4096, uint32] =
  var i = -1
  var counted:int = 0
  var skipped = newSeqOfCap[Record](skip_reads)
  for aln in bam:
    i += 1
    if i < skip_reads:
      skipped.add(aln.copy())
      continue
    else:
      skipped.setLen(0)
    if not aln.flag.proper_pair: continue
    if aln.isize < 0: continue
    if aln.isize > result.len: continue
    result[aln.isize].inc
    counted += 1
    if counted > n_reads: break

  if result.sum == 0:
    # mostly for debugging and testing on small bams.
    stderr.write_line "using first reads in fragment_length_distribution calculation as there were not enough"
    for aln in skipped:
      if not aln.flag.proper_pair: continue
      if aln.isize < 0 or aln.isize > result.len: continue
      result[aln.isize].inc

proc median(fragment_sizes: array[4096, uint32], pct:float=0.5): int =
  var n = sum(fragment_sizes)
  var count = 0'u32
  for i, cnt in fragment_sizes:
    count += cnt
    if count >=  uint32(0.5 + n.float / (1.0 / pct)):
      return i
  return fragment_sizes.len

proc inc[T](s:var Seq[T], enc:uint64) {.inline.} =
  s.A[enc.int].inc
  if s.imax == -1 or s.A[enc] > s.A[s.imax]:
    s.imax = enc.int

proc argmax[T](s: Seq[T]): uint64 {.inline.} =
  return s.imax.uint64

proc clear[T](s: var Seq[T]) {.inline.} =
  if s.imax == -1: return
  zeroMem(s.A[0].addr, sizeof(T) * len(s.A))
  s.imax = -1

iterator slide_by*(s:string, k: int): uint64 {.inline.} =
  ## given a string (DNA seq) yield the minimum kmer on the forward strand
  var base: char
  for i in countup(0, s.high - k + 1, k):
    var f = s[i..<i+k].encode()
    var kmin = f
    for j in 0..<k:
      base = s[i + j]
      f.forward_add(base, k)
      kmin = min(kmin, f)
    yield kmin


proc count(read: var string, k: int, count: var Seq[uint8]): int {.inline.} =
  # count the repeats of length k in read and return the most frequent
  count.clear
  for enc in read.slide_by(k):
    count.inc(enc)
  if count.imax == -1: return 0
  return count.A[count.imax].int

# This is the bottleneck for run time at the moment
proc get_repeat(read: var string, counts: var Seqs[uint8], repeat_count: var int, opts:Options): array[6, char] =
  if read.count('N') > 20: return

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
      var s = newString(k)
      counts[k].argmax.decode(s)
      if read.count(s) > (read.len.float * opts.proportion_repeat / k.float).int:
        copyMem(result[0].addr, s[0].addr, k)
        repeat_count = count
    elif count < (read.len.float * 0.12 / k.float).int:
      # e.g. for a 5 mer repeat, we should see some 2, 3, 4 mer reps and we can
      # bail if we do not. this is an optimization to avoid counting when we
      # cant possibly see a repeat.
      break

proc get_repeat(aln:Record, counts:var Seqs[uint8], repeat_count: var int, align_length: var int, opts:Options): array[6, char] =
  # returns blank array if nothing passes.
  var read = ""
  aln.sequence(read)
  align_length = len(read)
  result = read.get_repeat(counts, repeat_count, opts)

proc tostring(t:tread, targets: seq[Target]): string =
  var chrom = if t.tid == -1: "unknown" else: targets[t.tid].name
  var rep: string
  for v in t.repeat:
    if v == 0.char: continue
    rep.add(v)
  result = &"""{chrom}	{t.position}	{rep}	{t.split}	{t.repeat_count}"""
  when defined(debug):
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
  if aln.cigar.len > 0:
    if aln.cigar[0].op == CigarOp.soft_clip:
      align_length -= aln.cigar[0].len
    if aln.cigar[aln.cigar.len-1].op == CigarOp.soft_clip:
      align_length -= aln.cigar[aln.cigar.len-1].len

  result = tread(tid:aln.tid.int32,
                 position: aln.start.uint32,
                 repeat: rep,
                 flag: aln.flag,
                 repeat_count: repeat_count.uint8,
                 align_length: align_length.uint8,
                 split: Soft.none,
                 mapping_quality: aln.mapping_quality)
  when defined(debug):
    result.qname = aln.qname

type Cache = object
  tbl: TableRef[string, tread]
  cache: seq[tread]

proc add_soft(cache:var Cache, aln:Record, counts: var Seqs[uint8], opts:Options, read_repeat:array[6, char]) =
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
    var position = if cig_index == 0: (aln.start).uint32 else: (aln.stop).uint32
    cache.cache.add(tread(tid:aln.tid.int32,
                  position: position,
                  flag: aln.flag,
                  repeat: repeat,
                  repeat_count: repeat_count.uint8,
                  align_length: soft_seq.len.uint8,
                  split: if cig_index == 0: Soft.left else: Soft.right,
                  mapping_quality: aln.mapping_quality
                  ))
    when defined(debug):
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

proc add(cache:var Cache, aln:Record, counts: var Seqs[uint8], opts:Options) =
  doAssert not (aln.flag.secondary or aln.flag.supplementary)

  # Check if you have both reads
  if aln.after_mate:
    var added = false
    var mate:tread
    if not cache.tbl.take(aln.qname, mate): return

    # get mates together and see what, if anything should be added to the final
    # cache.

    var self = aln.to_tread(counts, opts)
    cache.add_soft(aln, counts, opts, self.repeat)
    if mate.repeat_count == 0'u8 and self.repeat_count == 0: return


    if mate.repeat_count > 0'u8:
      # mate is STR, mate is mapped well
      if mate.mapping_quality >= opts.min_mapq or mate.flag.proper_pair:
        added = true
        mate.position += uint32(mate.align_length.float / 2'f + 0.5) # Record position as middle of mate
        cache.cache.add(mate)
      # mate is STR, mate is mapped poorly, self is mapped well
      elif self.mapping_quality >= opts.min_mapq and not mate.flag.proper_pair:
        # self is right of mate, so the position subtracts th fragment length
        # Note fragment size is the external distance
        if self.flag.reverse:
          #   mate                  self
          #   =======>             <===========
          #   000000000000000000000000000000000 fragment length
          mate.position = self.position - opts.median_fragment_length.uint32 + self.align_length + uint32(mate.align_length.float / 2'f + 0.5)
        else:
          #   self                 mate
          #   =======>             <===========
          #   000000000000000000000000000000000 fragment length
          mate.position = self.position + opts.median_fragment_length.uint32 - uint32(mate.align_length.float / 2'f + 0.5)

        mate.tid = self.tid
        if mate.flag.should_reverse:
          mate.repeat.min_rev_complement
        added = true
        cache.cache.add(mate)

    if self.repeat_count > 0'u8:
      # self is STR, self is mapped well
      if self.mapping_quality >= opts.min_mapq or self.flag.proper_pair:
        added = true
        self.position += uint32(self.align_length.float / 2'f + 0.5)
        cache.cache.add(self)
      # self is STR, self is mapped poorly, mate is mapped well
      elif mate.mapping_quality >= opts.min_mapq and not self.flag.proper_pair:
        # self is right of mate, so the position subtracts the fragment length
        if mate.flag.reverse:
          #   self                 mate
          #   =======>             <===========
          #   000000000000000000000000000000000 fragment length
          self.position = mate.position + mate.align_length - opts.median_fragment_length.uint32 + uint32(mate.align_length.float / 2'f + 0.5)

        else:
          #   mate                  self
          #   =======>             <===========
          #   000000000000000000000000000000000 fragment length
          self.position = mate.position + opts.median_fragment_length.uint32 - uint32(self.align_length.float / 2'f + 0.5)
        self.tid = mate.tid
        if self.flag.should_reverse:
          self.repeat.min_rev_complement
        added = true
        cache.cache.add(self)

  else:
    var tr = aln.to_tread(counts, opts)
    cache.add_soft(aln, counts, opts, tr.repeat)
    doAssert not cache.tbl.hasKeyOrPut(aln.qname, tr), "error with read:" & aln.qname & " already in table as:" & $cache.tbl[aln.qname]

when isMainModule:
  import math




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
  var ibam:Bam
  var proportion_repeat = parseFloat(args.proportion_repeat)
  var skip_reads = parseInt(args.skip)

  if not open(ibam, args.bam, fai=args.fasta, threads=2):
    quit "couldn't open bam"

  var cram_opts = 8191 - SAM_RNAME.int - SAM_RGAUX.int - SAM_QUAL.int - SAM_SEQ.int
  discard ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, cram_opts)

  var frag_dist = ibam.fragment_length_distribution(skip_reads=skip_reads)
  var frag_median = frag_dist.median
  if args.verbose:
    stderr.write_line "Calculated median fragment length:", frag_median

  ibam.close()
  if not open(ibam, args.bam, fai=args.fasta, threads=2, index=true):
    quit "couldn't open bam"
  cram_opts = 8191 - SAM_RNAME.int - SAM_RGAUX.int - SAM_QUAL.int
  discard ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, cram_opts)


  var decodeds = newSeq[string](7)
  for i, s in decodeds.mpairs:
    decodeds[i] = newString(i)
  shallow(decodeds)

  var cache = Cache(tbl:newTable[string, tread](8192), cache: newSeqOfCap[tread](65556))
  var opts = Options(median_fragment_length: frag_median, proportion_repeat: proportion_repeat, min_mapq: 20'u8)

  var nreads = 0
  var counts = init[uint8]()
  for aln in ibam: #.query("14:92537254-92537477"):
    if aln.flag.secondary or aln.flag.supplementary: continue

    nreads.inc

    if nreads mod 10_000_000 == 0:
      var nrps = nreads.float64 / (cpuTime() - t0)
      if args.verbose:
        stderr.write_line $nreads, &" @ {aln.chrom}:{aln.start} {nrps:.1f} reads/sec tbl len: {cache.tbl.len} cache len: {cache.cache.len}"

    cache.add(aln, counts, opts)

  var reads_fh:File
  var bounds_fh:File
  if not open(reads_fh, args.output_prefix & "-reads.txt", mode=fmWrite):
    quit "couldn't open output file"
  if not open(bounds_fh, args.output_prefix & "-bounds.txt", mode=fmWrite):
    quit "couldn't open output file"

  reads_fh.write_line "chrom\tpos\tstr\tsoft_clip\tstr_count\tqname\tcluster_id" # print header
  var targets = ibam.hdr.targets
  var ci = 0
  for c in cache.cache.cluster(max_dist=frag_dist.median(0.98).uint32, min_supporting_reads=1):
    bounds_fh.write_line c.bounds.tostring(targets)
    for s in c.reads:
      reads_fh.write_line s.tostring(targets) & "\t" & $ci
    ci += 1

  reads_fh.close
  bounds_fh.close
  if args.verbose:
    stderr.write_line cache.tbl.len, " left in table"
    stderr.write_line &"wrote bounds to {args.output_prefix}-bounds.txt"
    stderr.write_line &"wrote reads to {args.output_prefix}-reads.txt"

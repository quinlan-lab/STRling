import kmer
import strutils
import hts/bam
import math
import tables
import strformat
import algorithm

{.push checks:off optimization:speed.}
iterator slide_by*(s:string, k: int): uint64 {.inline.} =
  ## given a string (DNA seq) yield the minimum kmer on the forward strand
  if k <= s.len:
    var base: char
    var f = s[0..<k].encode
    var kmin = f
    # note that we are just rotating the kmer here, not adding new bases
    for j in 0..<k:
        base = s[j]
        f.forward_add(base, k)
        kmin = min(f, kmin)
    yield kmin

    # after the first k, then we can use the forward add of k bases
    # to get to the next encode
    for i in countup(k, s.high - k + 1, k):
      for m in 0..<k:
        f.forward_add(s[i + m], k)
      kmin = f
      # then rotate the kmer
      for j in 0..<k:
          base = s[i + j]
          f.forward_add(base, k)
          kmin = min(f, kmin)
      yield kmin
{.pop.}

proc complement*(s:char): char {.inline.} =
    if s == 'C':
        return 'G'
    elif s == 'G':
        return 'C'
    elif s == 'A':
        return 'T'
    elif s == 'T':
        return 'A'
    else:
        return s

proc reverse_complement*(xs: string): string =
  result = newString(xs.len)
  for i, x in xs:
    # high == len - 1
    result[xs.high-i] = complement(x)

proc complement*(xs: string): string =
  result = newString(xs.len)
  for i, x in xs:
    # high == len - 1
    result[i] = complement(x)

proc min_rev_complement*(repeat: var array[6, char]) {.inline.} =
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
  s = s & s
  var mv = uint64(0) - 1'u64
  for m in s.slide_by(l):
    if m < mv:
      mv = m
  var ms = newString(l)
  mv.decode(ms)
  for i in 0..<l:
    repeat[i] = ms[i]

  # Return either the original or reverse complement as a string, 
  # whichever is first alphabetically


proc fragment_length_distribution*(bam:Bam, n_reads:int=2_000_000, skip_reads:int=100_000): array[4096, uint32] =
  var i = -1
  var counted:int = 0
  var skipped = newSeqOfCap[Record](skip_reads)
  for aln in bam:
    i += 1
    if not aln.flag.proper_pair: continue
    if aln.flag.supplementary or aln.flag.secondary: continue
    if aln.isize < 0: continue
    if aln.isize > result.len: continue
    if i < skip_reads:
      skipped.add(aln.copy())
      continue
    else:
      skipped.setLen(0)
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

type Seq*[T] = object
  imax*: int
  A*: seq[T]

type Seqs*[T] = array[7, Seq[T]]

type Options* = object
  median_fragment_length*: int
  proportion_repeat*: float
  min_mapq*: uint8
  min_support*: int
  window*: int
  targets*: seq[Target]

proc percentile*(fragment_sizes: array[4096, uint32], fragment_length:int): float =
  var total = sum(fragment_sizes)

  var s = 0
  for i, cnt in fragment_sizes:
    s += cnt.int
    if i >= fragment_length: break

  return s.float / total.float

proc median*(fragment_sizes: array[4096, uint32], pct:float=0.5): int =
  var n = sum(fragment_sizes)
  var count = 0'u32
  for i, cnt in fragment_sizes:
    count += cnt
    if count >=  uint32(0.5 + n.float / (1.0 / pct)):
      return i
  return fragment_sizes.len

proc median_depth*(D: seq[int]): int =
  var H = newSeq[int](1048)
  for d in D:
    when defined(debug):
      doAssert d >= 0
    H[min(d, H.high)] += 1 # Values greater than 1047 will be set to 1047 to avoid overflow
  var s: int
  for i, h in H:
    s += h
    if float(s) > float(D.len)/2.0:
      return i

proc mode*[T](xs: openArray[T]): T =
  var count = initCountTable[T]()
  return xs.toCountTable.largest.key

# Report the n most frequent keys in a CountTable, in decending order
proc most_frequent*[T](table: var CountTable[T], n: int): seq[T] =
  table.sort()
  if n > len(table):
    raise newException(IndexError, &"Insufficient keys in CountTable ({len(table)}) to report {n}")
  result = newseq[T](n)
  var i = 0
  for key, count in table:
    if i < n:
      result[i] = key
      inc i
    else:
      break

proc isNaN(v:float): bool {.inline.} =
  return v.classify == fcNaN

proc init*[T](): Seqs[T] =
  result = [
     Seq[T](A: newSeq[T](0)),
     Seq[T](A: newSeq[T](0)),
     Seq[T](A: newSeq[T](16)),
     Seq[T](A: newSeq[T](64)),
     Seq[T](A: newSeq[T](256)),
     Seq[T](A: newSeq[T](1024)),
     Seq[T](A: newSeq[T](4096)),
  ]

proc inc*[T](s:var Seq[T], enc:uint64) {.inline.} =
  s.A[enc.int].inc
  if s.imax == -1 or s.A[enc] > s.A[s.imax]:
    s.imax = enc.int

proc argmax*[T](s: Seq[T]): uint64 {.inline.} =
  return s.imax.uint64

proc clear*[T](s: var Seq[T]) {.inline.} =
  if s.imax == -1: return
  zeroMem(s.A[0].addr, sizeof(T) * len(s.A))
  s.imax = -1

proc count*(read: var string, k: int, count: var Seq[uint8]): int {.inline.} =
  # count the repeats of length k in read and return the most frequent
  count.clear
  for enc in read.slide_by(k):
    count.inc(enc)
  if count.imax == -1: return 0
  return count.A[count.imax].int

# Get a tid for a given chromosome name
proc get_tid*(name:string, targets: seq[Target]): int =
  for t in targets:
    if t.name == name:
      return t.tid
  return -1

proc reduce_repeat*(rep: var array[6, char]): int =
  ## change e.g. AA to A. result is the mutliplier, so AA to A becomes 2.
  ## CCC to C is 1.
  result = 1
  if rep[0] == '\0': return
  let seen = rep[0]
  for i in 1..<rep.len:
    if rep[i] == '\0': break
    if rep[i] != seen: return
  # if we got here, we have a homopolymer
  for i in 1..<rep.len:
    if rep[i] == '\0': break
    result.inc
    rep[i] = '\0'

# This is the bottleneck for run time at the moment
proc get_repeat*(read: var string, counts: var Seqs[uint8], repeat_count: var int, opts:Options, debug:bool=false): array[6, char] =
  repeat_count = 0
  if read.count('N') > 20: return
  var s = newString(6)

  var best_score: int = -1
  for k in 2..6:
    var count = read.count(k, counts[k])
    s = newString(k)
    counts[k].argmax.decode(s)
    var score = count * k
    if debug:
      echo k, " ", score, " ", s, " actual count:", read.count(s), " est count:", count
    if score <= best_score:
      if count < (read.len.float * 0.12 / k.float).int:
        break
      continue
    count = read.count(s)
    score = count * k
    if score < best_score: continue

    best_score = score
    if count > (read.len.float * opts.proportion_repeat / k.float).int:
      # now check the actual string because the kmer method can't track phase
      if count >= (read.len.float * opts.proportion_repeat / k.float).int:
        copyMem(result.addr, s[0].addr, k)
        repeat_count = count
        if repeat_count > 0 and result[0] == '\0':
          quit "bad:" & $k & " " & $result & " " & "kmer:" & s
    elif count < (read.len.float * 0.12 / k.float).int:
      # e.g. for a 5 mer repeat, we should see some 2, 3, 4 mer reps and we can
      # bail if we do not. this is an optimization to avoid counting when we
      # cant possibly see a repeat.
      break

  repeat_count *= reduce_repeat(result)

proc tostring*(a:array[6, char]): string =
  for c in a:
    if c == 0.char: return
    result.add(c)

# Represent an STR repeat unit as an array
proc as_array*(s:string): array[6, char] =
  when defined(debug):
    doAssert s.len <= 6
    for x in s:
      doAssert x in "ATCG"
  for i, c in s:
    result[i] = c

# Get a chromosome name for a given tid
proc get_chrom*(tid:int, targets: seq[Target]): string =
  for t in targets:
    if t.tid == tid:
      return t.name

proc `<`(a: array[6, char], b: array[6, char]): bool {.inline.} =
  if a[0] != b[0]:
    return a[0] < b[0]
  if a[1] != b[1]:
    return a[1] < b[1]
  if a[2] != b[2]:
    return a[2] < b[2]
  if a[3] != b[3]:
    return a[3] < b[3]
  if a[4] != b[4]:
    return a[4] < b[4]
  return a[5] < b[5]

proc canonical_repeat*(repeat: array[6, char]): array[6, char] =
  for i in 0..<6:
    result[i] = repeat[i]
  result.min_rev_complement
  if result < repeat:
    return
  return repeat

proc canonical_repeat*(repeat: string): string =
  var forward = ['\0', '\0', '\0', '\0', '\0', '\0']
  for i in 0..<len(repeat):
    forward[i] = repeat[i]
  result = forward.canonical_repeat.tostring


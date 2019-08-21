import kmer
import math
import hts/bam

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


type Seq*[T] = object
  imax*: int
  A*: seq[T]

type Seqs*[T] = array[7, Seq[T]]

type Options* = object
  median_fragment_length*: int
  proportion_repeat*: float
  min_mapq*: uint8
  min_support*: int

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



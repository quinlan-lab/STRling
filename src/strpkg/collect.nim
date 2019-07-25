import hts/bam
import ./cluster

type Support* = object
  SpanningFragmentLength*: uint32
  SpanningReadRepeatCount*: uint8
  SpanningReadCigarInsertionCount*: uint8
  SpanningReadCigarDeletionCount*: uint8


proc spanning_fragment*(L:Record, R:Record, bounds:Bounds, support:var Support): bool =
  doAssert L.start <= R.start
  if L.start < bounds.left.int and R.stop > bounds.right.int:
    support.SpanningFragmentLength = L.isize.abs.uint32
    result = true

proc spanning_read(A:Record, bounds:Bounds, support: var Support): bool =

  if A.start < bounds.left.int and A.stop > bounds.right.int:

    for cig in A.cigar:
      if cig.op == Cigarop.insert:
        support.SpanningReadCigarInsertionCount += 1
      if cig.op == Cigarop.deletion:
        support.SpanningReadCigarDeletionCount += 1
    result = true


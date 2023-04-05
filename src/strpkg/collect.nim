import hts/bam
import strutils
import algorithm
import strformat
import tables
import math
import ./cluster
import ./utils

type SupportType* {.pure.} = enum
    SpanningFragment # A spanning pair of reads definine a fragment
    SpanningRead # A read that spans the locus
    OverlappingRead # A read that overlaps but does not span the locus (e.g. is soft-clipped)

type Support* = object
  Type*: SupportType
  
  # if SpanningFragmentLength is non-zero, we are looking at a Fragment
  # and all other fields will be 0
  SpanningFragmentLength*: uint32
  SpanningFragmentPercentile*: float

  # the *Read* fields are for non-fragments though they may be 0 if
  # for example, the reference has no STR and we're looking at the reference
  # allele.
  SpanningReadRepeatCount*: uint8
  # counting number of deletion ops.
  SpanningReadCigarInsertionLen*: uint8
  SpanningReadCigarDeletionLen*: uint8
  repeat*: string
  qname*: string

proc tostring*(s:Support, b:Bounds, chrom:string): string =
  result = &"{chrom}\t{b.left}\t{b.right}\t{s.Type}\t{s.SpanningFragmentLength}\t{s.SpanningFragmentPercentile}\t{s.SpanningReadRepeatCount}\t{s.SpanningReadCigarInsertionLen}\t{s.SpanningReadCigarDeletionLen}\t{s.repeat}\t{s.qname}"

proc spanning_fragment*(L:Record, R:Record, bounds:Bounds, support:var Support, frag_sizes: array[4096, uint32]): bool =
  doAssert L.start <= R.start, $(L, R, bounds)
  var bound_width = bounds.right.int - bounds.left.int
  var slop = bounds.repeat.len - 1
  if bound_width < 5: # Add extra slop for very small bounds
    slop += (5 - bound_width)
  if L.start < (bounds.left.int - slop) and R.stop > (bounds.right.int + slop):
    support.Type = SupportType.SpanningFragment
    support.SpanningFragmentLength = max(1'u32, L.isize.abs.uint32)
    support.SpanningFragmentPercentile = frag_sizes.percentile(support.SpanningFragmentLength.int)
    support.repeat = bounds.repeat
    support.qname = L.qname
    result = true

proc find_read_position(A: Record, position:int): int =
  var
    r_off = A.start.int
    q_off = 0
  result = -1

  for op in A.cigar:
    if r_off > position: return -1
    let cons = op.consumes
    if cons.query:
      q_off += op.len
    if cons.reference:
      r_off += op.len
    if r_off < position: continue

    var over = r_off - position
    if over > q_off: return -1

    if not cons.query:
      return -1

    return q_off - over


proc count(A: Record, bounds:Bounds): int =
  ## given a read, count the repeat units in that read that are in the
  ## reference offset provided in bounds.
  var dna:string
  if bounds.right < bounds.left: return 0
  A.sequence(dna)

  var read_left = A.find_read_position(bounds.left.int)
  var read_right = A.find_read_position(bounds.right.int)
  if read_left >= 0 and read_right < 0:
    read_right = dna.len
  if read_left < 0 and read_right < 0:
    return 0

  if read_left < 0: read_left = 0
  let S = dna[read_left..<read_right]
  result = S.count(bounds.repeat)
  if result < int(S.len.float * 0.7 / bounds.repeat.len.float):
    result = 0

# Returns true if a read overlaps the bounds. Also records if the read spans
# and counts STRs
proc overlapping_read*(A:Record, bounds:Bounds, support: var Support): bool =
  var bound_width = bounds.right.int - bounds.left.int
  var slop = bounds.repeat.len - 1
  if bound_width < 5: # Add extra slop for very small bounds
    slop += (5 - bound_width)
  # Check if read overlaps bounds
  if A.overlaps(bounds):
    result = true
    support.Type = SupportType.OverlappingRead
    # Count STRs
    support.SpanningReadRepeatCount = A.count(bounds).uint8
    support.qname = A.qname

    # Check if read completely spans bounds
    if A.start < (bounds.left.int - slop) and A.stop > (bounds.right.int + slop):
      support.Type = SupportType.SpanningRead
      for cig in A.cigar:
        if cig.op == Cigarop.insert:
          support.SpanningReadCigarInsertionLen += cig.len.uint8
        if cig.op == Cigarop.deletion:
          support.SpanningReadCigarDeletionLen += cig.len.uint8

proc estimate_size*(spanners: seq[Support], frag_sizes: array[4096, uint32]): int =
  var small_sizes = newSeq[uint32]()
  for s in spanners:
    if s.SpanningFragmentLength > 0'u32 and s.SpanningFragmentPercentile < 0.01:
      small_sizes.add(s.SpanningFragmentLength)
  if small_sizes.len == 0: return -1
  sort(small_sizes)
  var s = small_sizes[int(small_sizes.high/2)]
  return frag_sizes.median - s.int

import ./spanning

proc spanners*(b:Bam, bounds:Bounds, window:int, frag_sizes: array[4096, uint32], min_mapq:uint8=20, max_size:int=5000): tuple[support: seq[Support], median_depth: int, expected_spanners: float32] =
  var pairs = newTable[string, seq[Record]]()
  doAssert left <= right
  var window_left = bounds.left.int - window
  var window_right = bounds.right.int + window
  var cd = frag_sizes.cumulative
  var depths = newSeq[int](window_right - window_left) # depth per base within the bounds

  var expected_spanners_by_qname = initTable[string, float]()
  for aln in b.query(bounds.tid.int, max(0, window_left), window_right):
     if aln.flag.secondary or aln.flag.supplementary or aln.flag.dup: continue
     if aln.mapping_quality < min_mapq: continue
     var prob = cd.expected_spanning_probability(aln, bounds.left.int, bounds.right.int)
     if prob > 0:

       if aln.qname in expected_spanners_by_qname:
         # if already present, we take average of estimate from each pair
         expected_spanners_by_qname[aln.qname] = 0.5 * (expected_spanners_by_qname[aln.qname] + cd.expected_spanning_probability(aln, bounds.left.int, bounds.right.int))
       else:
         expected_spanners_by_qname[aln.qname] = cd.expected_spanning_probability(aln, bounds.left.int, bounds.right.int)

     depths[max(0, aln.start - window_left - 1)] += 1
     depths[min(depths.high, aln.stop - window_left - 1)] -= 1

     var s:Support
     if aln.overlapping_read(bounds, s):
       result.support.add(s)
     if aln.tid != aln.mate_tid: continue
     if aln.isize.abs > max_size: continue

     # we could check here if we have 2 in the pair, process immediately, and remove.
     # but deleting from a table is not cheap so current trade-off memory for speed.
     # should instead add a reduced object with chrom, sequence, cigar
     pairs.mgetOrPut(aln.qname, @[]).add(aln.copy())
     when defined(debug):
      if pairs.len mod 5000 == 0:
        stderr.write_line "memory usage. number of pairs:", pairs.len, " bounds:", bounds
     if pairs.len > 20_000:
       when defined(debug):
         stderr.write_line "high-depth for bounds:", bounds, " skipping"
       return (@[], -1, 0'f32)

  for v in expected_spanners_by_qname.values:
    result.expected_spanners += v

  for qname, pair in pairs:
    if len(pair) != 2: continue
    var s: Support
    if spanning_fragment(pair[0], pair[1], bounds, s, frag_sizes):
      result.support.add(s)

  depths.cumsum
  result.median_depth = median_depth(depths)


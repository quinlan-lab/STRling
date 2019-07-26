import hts/bam
import strutils
import algorithm
import strformat
import tables
import ./cluster
import ./utils

type Support* = object
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
  when defined(debug) or defined(qname):
    qname: string

proc tostring*(s:Support, b:Bounds, chrom:string): string =
  result = &"{chrom}\t{b.left}\t{b.right}\t{s.SpanningFragmentLength}\t{s.SpanningFragmentPercentile}\t{s.SpanningReadRepeatCount}\t{s.SpanningReadCigarInsertionLen}\t{s.SpanningReadCigarDeletionLen}"
  when defined(debug) or defined(qname):
    result &= "\t" & s.qname

proc spanning_fragment*(L:Record, R:Record, bounds:Bounds, support:var Support, frag_sizes: array[4096, uint32]): bool =
  doAssert L.start <= R.start
  if L.start < bounds.left.int and R.stop > bounds.right.int:
    support.SpanningFragmentLength = max(1'u32, L.isize.abs.uint32)
    support.SpanningFragmentPercentile = frag_sizes.percentile(support.SpanningFragmentLength.int)
    when defined(debug) or defined(qname):
      support.qname = L.qname
    result = true

proc find_read_position(A: Record, position:int): int =
  var
    r_off = A.start
    q_off = 0

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

    return q_off - over


proc count(A: Record, bounds:Bounds): int =
  ## given a read, count the repeat units in that read that are in the
  ## reference offset provided in bounds.
  var dna:string
  A.sequence(dna)

  var read_left = A.find_read_position(bounds.left.int)
  var read_right = A.find_read_position(bounds.right.int)
  if read_left < 0: read_left = 0
  if read_right < 0: read_right = 1

  return dna[read_left..read_right].count(bounds.repeat)

proc spanning_read*(A:Record, bounds:Bounds, support: var Support): bool =

  if A.start < bounds.left.int and A.stop > bounds.right.int:
    # just do a count directly
    support.SpanningReadRepeatCount = A.count(bounds).uint8
    when defined(debug) or defined(qname):
      support.qname = A.qname

    for cig in A.cigar:
      # TODO: check if the sequence of this cigar op matches the repeat unit
      if cig.op == Cigarop.insert:
        support.SpanningReadCigarInsertionLen += cig.len.uint8
      if cig.op == Cigarop.deletion:
        support.SpanningReadCigarDeletionLen += cig.len.uint8
    result = true

proc estimate_size*(spanners: seq[Support], frag_sizes: array[4096, uint32]): int =
  var small_sizes = newSeq[uint32]()
  for s in spanners:
    if s.SpanningFragmentLength > 0'u32 and s.SpanningFragmentPercentile < 0.01:
      small_sizes.add(s.SpanningFragmentLength)
  if small_sizes.len == 0: return -1
  sort(small_sizes)
  var s = small_sizes[int(small_sizes.high/2)]
  return frag_sizes.median - s.int

proc spanners*(b:Bam, bounds:Bounds, window:int, frag_sizes: array[4096, uint32], min_mapq:uint8=20): seq[Support] =
  var pairs = newTable[string, seq[Record]]()
  for aln in b.query(bounds.tid.int, max(0, bounds.left.int - window), bounds.right.int + window):
     if aln.flag.secondary or aln.flag.supplementary: continue
     if aln.mapping_quality < min_mapq: continue

     var s:Support
     if aln.spanning_read(bounds, s):
       result.add(s)

     pairs.mgetOrPut(aln.qname, @[]).add(aln.copy())

  for qname, pair in pairs:
    if len(pair) != 2: continue
    var s: Support
    if spanning_fragment(pair[0], pair[1], bounds, s, frag_sizes):
      result.add(s)

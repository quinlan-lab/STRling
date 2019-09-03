import algorithm
import math
import strformat
import tables
import itertools
import hts/bam
import strutils
import utils
import msgpack4nim
import msgpack4collection

type Soft* {.size:1, pure.} = enum
  left
  right
  both
  none

# Data structure storing information about each read that looks like an STR
type tread* = object
  tid*: int32
  position*: uint32
  repeat*: array[6, char]
  flag*: Flag
  split*: Soft
  mapping_quality*: uint8
  repeat_count*: uint8
  align_length*: uint8
  when defined(debug) or defined(qname):
    qname*: string

proc pack_type*[ByteStream](s: ByteStream, x: tread) =
  s.pack(x.tid)
  s.pack(x.position)
  s.pack(x.repeat)
  s.pack(x.flag.uint16)
  s.pack(x.split.uint8)
  s.pack(x.mapping_quality)
  s.pack(x.repeat_count)
  s.pack(x.align_length)
  var L:uint32 = 0
  when defined(debug) or defined(qname):
    L = x.qname.len.uint32
    s.pack(L)
    s.pack(x.qname)
  else:
    s.pack(L)

proc unpack_type*[ByteStream](s: ByteStream, x: var tread) =
  s.unpack(x.tid)
  s.unpack(x.position)
  s.unpack(x.repeat)
  var f:uint16
  s.unpack(f)
  x.flag = Flag(f)
  var split:uint8
  s.unpack(split)
  x.split = Soft(split)
  s.unpack(x.mapping_quality)
  s.unpack(x.repeat_count)
  s.unpack(x.align_length)
  var L:uint32 = 0
  var qname: string
  s.unpack(L)
  qname = newString(L)
  s.unpack(qname)
  when defined(debug) or defined(qname):
    x.qname = qname


type Cluster* = object
  reads*: seq[tread]
  left: int
  right:int

const mediani = 9

proc posmed(cl:Cluster, n:int=5): uint32 =
  ## posmed is the median of the first n positions in the cluster
  let mid = int(min(mediani, cl.reads.len) / 2 - 0.5)
  return cl.reads[mid].position

proc bytidrep(t:tread): tuple[tid:int32, repeat:array[6, char]] =
  return (t.tid, t.repeat)

# Sorts the reads by chromosome (tid) then repeat unit, then by position
proc tread_cmp(a: tread, b:tread): int =
  if a.tid != b.tid: return cmp(a.tid, b.tid)
  for i in 0..<6:
    if a.repeat[i] != b.repeat[i]:
      return cmp(a.repeat[i], b.repeat[i])
  return cmp(a.position, b.position)

type Bounds* = object
  tid*: int32
  left*: uint32
  right*: uint32
  center_mass*: uint32
  n_left*: uint16
  n_right*: uint16
  n_total*: uint16
  repeat*: string
  name*: string
  left_exact*: bool
  right_exact*:bool


proc `==`*(a,b: Bounds): bool =
  if (a.tid == b.tid) and (a.left == b.left) and (a.right == b.right) and (a.repeat == b.repeat):
    return true

# Check if two Bounds overlap. Assumes left <= right in both Bounds
proc overlaps*(a,b: Bounds): bool =
  if a.tid == b.tid and a.repeat == b.repeat:
    var ileft = max(a.left, b.left) #lower bound of intersection interval
    var iright = min(a.right, b.right) #upper bound of intersection interval
    return ileft <= iright #interval non-empty?

# Parse single line of a an STR loci file
proc parse_bounds*(l:string, targets: seq[Target]): Bounds =
  var l_split = l.splitWhitespace()
  if len(l_split) == 4:
    discard
  elif len(l_split) == 5:
    result.name = l_split[4]
  else:
    quit fmt"Error reading loci bed file. Expected 4 or 5 fields and got {len(l_split)} on line: {l}"
  result.tid = int32(get_tid(l_split[0], targets))
  result.left = uint32(parseInt(l_split[1]))
  result.right = uint32(parseInt(l_split[2]))
  result.repeat = l_split[3]
  doAssert(result.left <= result.right)

# Parse an STR loci bed file
proc parse_loci*(f:string, targets: seq[Target]): seq[Bounds] =
  for line in lines f:
    result.add(parse_bounds(string(line), targets))

# Find the bounds of the STR in the reference genome
proc bounds*(cl:Cluster): Bounds =

  var lefts = initCountTable[uint32](8)
  var rights = initCountTable[uint32](8)
  result.left = cl.left.uint32
  result.right = cl.right.uint32

  var posns = newSeqOfCap[uint32](cl.reads.len)
  for c in cl.reads[0].repeat:
    if c == 0.char: break
    result.repeat.add(c)
  result.tid = cl.reads[0].tid
  doAssert cl.reads.len <= uint16.high.int, ("got too many reads for cluster with first read:" & $cl.reads[0])

  for r in cl.reads:
    result.n_total.inc
    if r.split == Soft.right:
      # soft-clip of left indicates right end of variant
      rights.inc(r.position)
      result.n_right.inc
    elif r.split == Soft.left:
      lefts.inc(r.position)
      result.n_left.inc
    else:
      posns.add(r.position)
  if posns.len > 0:
    result.center_mass = posns[int(posns.len / 2)]
  # require at least 50% of lefts to have the same position.
  if lefts.len > 0 and lefts.largest.val.float > 0.5 * lefts.len.float:
    var ll = lefts.largest
    result.left = ll.key
    # only set to exact if most frequent is at least 2. and above, we've also
    # required to have > 50% of left-splits at this exact location.
    result.left_exact = ll.val > 1
  else:
    result.left = result.center_mass
  if rights.len > 0 and rights.largest.val.float > 0.5 * rights.len.float:
    var rr = rights.largest
    result.right = rr.key
    result.right_exact = rr.val > 1
  else:
    result.right = result.center_mass

  # if we have soft clips in left, but not right AND the bounds are
  # flipped set right == left + 1. (and vice-versa for right <- left)
  if result.right < result.left:
    if (not result.right_exact) and result.left_exact:
      result.right = result.left + 1
    elif (not result.left_exact) and result.right_exact:
      result.left = result.right - 1

  if result.left > result.right:
    result.left = cl.left.uint32
    result.right = cl.right.uint32

proc trim(cl:var Cluster, max_dist:uint32) =
  if cl.reads.len == 0: return
  # drop stuff from start of cluster that is now outside the expected distance
  var lo = max(0, cl.posmed(mediani).int - max_dist.int).uint32
  while len(cl.reads) > 1 and cl.reads[0].position < lo:
    cl.reads = cl.reads[1..cl.reads.high]

proc tostring*(b:Bounds, targets: seq[Target]): string =
  return &"{targets[b.tid].name}\t{b.left}\t{b.right}\t{b.center_mass}\t{b.n_left}\t{b.n_right}\t{b.n_total}\t{b.repeat}\t{b.name}"

proc tostring*(c:Cluster, targets: seq[Target]): string =
  var rep: string
  for v in c.reads[0].repeat:
    if v == 0.char: continue
    rep.add(v)
  return &"{targets[c.reads[0].tid].name}\t{c.reads[0].position}\t{c.reads[c.reads.high].position}\t{c.reads.len}\t{rep}"

proc has_anchor(reads: seq[tread]): bool =
  # we can get false clusters driven only soft-clips
  # and we can remove them by noting that there are no
  # anchoring reads (mates mapping well nearby)
  for r in reads:
    if r.split == Soft.none: return true
  return false

type pair = object
  left: int
  right: int

proc fix_overlapping_bounds(strong_pairs: var seq[pair]) =
  var to_rm: seq[int]
  for i in 1..strong_pairs.high:
    var a = strong_pairs[i - 1]
    var b = strong_pairs[i]
    #echo "before:", a, " ", b
    #doAssert a.left <= b.left, &"{a} should be left of {b} at {i}"
    if a.right <= b.left: continue
    var mid = int(0.5 + (a.right.float + b.left.float) / 2.0)
    strong_pairs[i - 1].right = mid
    var R = strong_pairs[i - 1]
    if R.right < R.left:
      strong_pairs[i - 1].left = R.right - 1

    var L = strong_pairs[i]
    strong_pairs[i].left = mid
    if L.right < L.left:
      strong_pairs[i].right = L.left + 1

    if strong_pairs[i-1].right > strong_pairs[i].left:
      to_rm.add(i-1)
      strong_pairs[i] = b

  for i in to_rm.reversed:
    strong_pairs.delete(i)

proc count_lefts_and_rights(tandems: seq[tread]): tuple[lefts: CountTableRef[uint32], rights: CountTableRef[uint32]] =
  result.lefts = newCountTable[uint32]()
  result.rights = newCountTable[uint32]()
  for i in 0..tandems.high:
    var t = tandems[i]
    if t.split == Soft.left:
      result.lefts.inc(t.position)
    elif t.split == Soft.right:
      result.rights.inc(t.position)

proc find_strong(sites: CountTableRef[uint32], strong_soft_cutoff:int, min_dist:uint32=3): seq[uint32] =
  # find sites with at least strong_soft_cutoff clips at the same position and
  # require each position to be at least min_dist away from a previous one.
  result = newSeqOfCap[uint32](2048)

  for pos, cnt in sites:
    if cnt >= strong_soft_cutoff:
      if result.len > 0 and pos - result[result.high] < min_dist:
        continue
      else:
        result.add(pos)
  sort(result)

proc cluster(strong_pairs:var seq[pair], tandems: var seq[tread]): seq[Cluster] =
  # given bounds, generate clusters
  for p in strong_pairs:
    if tandems.len == 0: return
    var ilo = tandems.lowerBound(p.left, proc(a: tread, i: int): int =
      return cmp(a.position.int, p.left)
    )
    if ilo < 0: ilo = 0

    var ihi = tandems.upperBound(p.right, proc(a: tread, i: int): int =
      return cmp(a.position.int, p.right)
    )
    if ihi < 0: ihi = 0 #tandems.len
    try:
      var candidates = tandems[ilo..<ihi]
      result.add(Cluster(reads:candidates, left: p.left, right: p.right))
      tandems = tandems[0..<ilo] & tandems[ihi..tandems.high]
    except:
      echo &"ERROR with pair: {p}"
      echo &"ilo:{ilo} ihi:{ihi}"
      raise getCurrentException()

iterator gen_strong_single_side(tandems: var seq[tread], max_dist:uint32, min_supporting_reads:int, strong_soft_cutoff:int): Cluster =
  var x = tandems[0].tid
  for t in tandems:
    doAssert t.tid == x
  var (lefts, rights) = tandems.count_lefts_and_rights
  var strong_lefts = lefts.find_strong(strong_soft_cutoff)
  var strong_rights = rights.find_strong(strong_soft_cutoff)

  # we call these pairs, but the bounds are derived from just a single side.
  var pairs = newSeq[pair]()
  for p in strong_lefts:
    pairs.add(pair(left: max(0, p.int - max_dist.int), right: p.int + max_dist.int))
  for p in strong_rights:
    pairs.add(pair(left: max(0, p.int - max_dist.int), right: p.int + max_dist.int))
  pairs.sort do (a, b: pair) -> int:
    result = cmp(a.left.int, b.left.int)
    if result == 0:
      result = cmp(a.right.int, b.right.int)
  pairs.fix_overlapping_bounds
  for c in pairs.cluster(tandems):
    yield c

iterator gen_strong_pairs(tandems: var seq[tread], max_dist:uint32, min_supporting_reads:int, strong_soft_cutoff:int): Cluster =

  # count-tables with pos-> count
  var (lefts, rights) = tandems.count_lefts_and_rights

  # strong_left/rights are seq postion for each position with >=
  # strong_soft_cutoff
  var strong_lefts = lefts.find_strong(strong_soft_cutoff)
  var strong_rights = rights.find_strong(strong_soft_cutoff)

  var strong_pairs = newSeq[pair]()
  if strong_rights.len > 0:

    for left in strong_lefts:
      # check if we have a right nearby
      var right = strong_rights[min(strong_rights.high, strong_rights.upperBound(left))]
      # find a right within 100 bases. TODO: make this a(n internal) parameter
      if right.int - left.int > 100:
        continue

      if left > right:
        #echo &"ERROR: left:{left} > right:{right} . upperBounds was: {strong_rights.upperBound(left)}"
        continue

      strong_pairs.add(pair(left: max(0, left.int - max_dist.int), right: right.int + max_dist.int))

    # adjust so bounds are non-overlapping
  for p in strong_pairs:
    doAssert p.left < p.right, $p
  strong_pairs.fix_overlapping_bounds
  for p in strong_pairs:
    doAssert p.left < p.right, $p
  for c in strong_pairs.cluster(tandems):
    yield c

iterator cluster_single(tandems: var seq[tread], max_dist:uint32, min_supporting_reads:int=5, strong_soft_cutoff:int=4): Cluster =
  # tandems are a single repeat unit from a single chromosome.
  if tandems[0].tid < 0:
    stderr.write_line "yielding " & $tandems.len & " unplaced reads with repaeat: " & $tandems[0].repeat
    yield Cluster(reads: tandems)
    tandems.setLen(0)
  else:
    var k = 0
    var L = tandems.len
    for c in gen_strong_pairs(tandems, max_dist, min_supporting_reads, strong_soft_cutoff):
      k += c.reads.len
      yield c
    doAssert tandems.len == L - k
    for c in gen_strong_single_side(tandems, max_dist, min_supporting_reads, strong_soft_cutoff):
      k += c.reads.len
      yield c
    doAssert tandems.len == L - k

    # TODO: cluster those with no or little split support.

iterator cluster*(tandems: var seq[tread], max_dist:uint32, min_supporting_reads:int=5, strong_soft_cutoff:int=4): Cluster =
  tandems.sort(tread_cmp)

  for group in groupby(tandems, bytidrep):
    # reps are on same chromosome and have same repeat unit
    var reps: seq[tread] = group.v
    # first we try clustering based on the split reads.
    for cluster in reps.cluster_single(max_dist, min_supporting_reads): yield cluster

    # TODO: continue here?
    #[

    if reps[0].tid < 0:
      stderr.write_line "yielding " & $reps.len & " unplaced reads with repaeat: " & $reps[0].repeat
      yield Cluster(reads: reps)
      continue

    var i = 0
    var c:Cluster
    while i < reps.len:
      # start a new cluster
      var it = reps[i]
      c = Cluster(reads: @[it])
      i += 1 # increment i even if we dont enter the loop
      for j in i..reps.high:
        # add any tread that's close enough. we use 2 * max_dist because
        # we can have reads a fragment length away from the middle of the event
        # from either direction.
        if reps[j].position <= c.posmed(mediani) + 2'u32 * max_dist:
          c.reads.add(reps[j])
          i = j + 1
          continue

        # remove stuff (at start of cluster) that's now too far away.
        c.trim(2'u32 * max_dist)
        if c.reads.len >= min_supporting_reads and c.reads.has_anchor:
          yield c
          c = Cluster()
        # increment i to past last j and break out of this cluster
        break

    c.trim(max_dist)
    if c.reads.len >= min_supporting_reads and c.reads.has_anchor:
      yield c
    ]#

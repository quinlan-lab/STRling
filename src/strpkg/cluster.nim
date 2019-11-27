import algorithm
import msgpack4nim
import msgpack4collection
import math
import strformat
import tables
import itertools
import hts/bam
import strutils
import utils

type Soft* {.size:1, pure.} = enum
  left ## the left clipped portion of the read is repetitive
  right ## the right clipped portion of the read is repetitive
  both ## both sides of read are clipped and both are repetitive
  none
  none_right ## looking at main part of read which is soft-clipped on the right
  none_left ## looking at main part of read which is soft-clipped on the left

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
  L = x.qname.len.uint32
  s.pack(L)
  s.pack(x.qname)

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
  if L > 0'u32:
    s.unpack(qname)
  x.qname = qname


type Cluster* = object
  reads*: seq[tread]

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
  force_report*: bool

proc `==`(a,b: Bounds): bool =
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
  for x in result.repeat:
    if x notin "ATCG":
      quit fmt"Error reading loci bed file. Expected DNA (ATCG only) in the 4th field, and got an unexpected character on line: {l}"
  doAssert(result.left <= result.right)

# Parse an STR loci bed file
proc parse_loci*(f:string, targets: seq[Target]): seq[Bounds] =
  for line in lines f:
    result.add(parse_bounds(string(line), targets))

# Find the bounds of the STR in the reference genome
proc bounds*(cl:Cluster): Bounds =

  var lefts = initCountTable[uint32](8)
  var rights = initCountTable[uint32](8)

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
    posns.add(r.position)

  if lefts.len > 0:
    var ll = lefts.largest
    if ll.val > 1:
      result.left = ll.key
  if rights.len > 0:
    var rr = rights.largest
    if rr.val > 1:
      result.right = rr.key

  if posns.len > 0:
    result.center_mass = posns[int(posns.len / 2)]
    if result.left == 0:
      result.left = result.center_mass
    if result.right == 0:
      result.right = result.left + 1
  else:
    if result.right == 0:
      result.right = result.left + 1
    if result.left == 0:
      result.left = result.right - 1

  # finally, here we can have the situation where we set the
  # left based on center-mass and it's larger than right
  # so we just swap them
  if result.left >= result.right:
    if result.n_left > 0'u16 and result.n_right > 0'u16:
      swap(result.left, result.right)
    else:
      result.left = result.right - 1

proc trim(cl:var Cluster, max_dist:uint32) =
  if cl.reads.len == 0: return
  # drop stuff from start of cluster that is now outside the expected distance
  var lo = max(0, cl.posmed(mediani).int - max_dist.int).uint32
  while len(cl.reads) > 1 and cl.reads[0].position < lo:
    cl.reads = cl.reads[1..cl.reads.high]

proc tostring*(b:Bounds, targets: seq[Target]): string =
  return &"{targets[b.tid].name}\t{b.left}\t{b.right}\t{b.repeat}\t{b.name}\t{b.center_mass}\t{b.n_left}\t{b.n_right}\t{b.n_total}"

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

iterator split_cluster*(c:Cluster, min_supporting_reads:int): Cluster =
  ## results for splitting
  ## 1. given a group of rights with min_supporting reads that is also the
  ## largest peak of rights in the current cluster *and* there is a peak
  ## of lefts that follows, we split between them to make 2 clusters.
  var should_split = false
  var lefts = initCountTable[uint32](8)
  var rights = initCountTable[uint32](8)
  for r in c.reads:
    if r.split == Soft.left: lefts.inc(r.position)
    elif r.split == Soft.right: rights.inc(r.position)

  if rights.len == 0 or lefts.len == 0:
    yield c

  else:

    var rl = rights.largest
    var ll = lefts.largest
    # NOTE this currently looks only at the largest.
    if rl.key < ll.key and rl.val > min_supporting_reads and ll.val > min_supporting_reads and ll.val.float / lefts.len.float > 0.5 and rl.val.float / rights.len.float > 0.5:
      var mid = uint32(0.5 + (rl.key.float + ll.key.float) / 2.0)
      var c1 = Cluster()
      var c2 = Cluster()
      for r in c.reads:
        if r.position < mid:
          c1.reads.add(r)
        else:
          c2.reads.add(r)

      yield c1
      yield c2

    else:
      yield c


iterator trcluster*(reps: seq[tread], max_dist:uint32, min_supporting_reads:int): Cluster =
  var i = 0
  var c:Cluster
  while i < reps.len:
    # start a new cluster
    var it = reps[i]
    c = Cluster(reads: @[it])
    i += 1 # increment i even if we dont enter the loop
    for j in i..reps.high:
      # add any tread that's close enough.
      # we can have reads a fragment0length away from the middle of the event
      # from either direction. add 100 for max event length
      if reps[j].position <= c.posmed(mediani) + max_dist + 100:
        c.reads.add(reps[j])
        i = j + 1
        continue

      # remove stuff (at start of cluster) that's now too far away.
      c.trim(max_dist + 100)
      if c.reads.len >= min_supporting_reads and c.reads.has_anchor:
        for sc in c.split_cluster(min_supporting_reads): yield sc
        c = Cluster()
      # increment i to past last j and break out of this cluster
      break

  c.trim(max_dist + 100)
  if c.reads.len >= min_supporting_reads and c.reads.has_anchor:
    for sc in c.split_cluster(min_supporting_reads): yield sc


iterator cluster*(tandems: var seq[tread], max_dist:uint32, min_supporting_reads:int=5): Cluster =
  tandems.sort(tread_cmp)

  for group in groupby(tandems, bytidrep):
    # reps are on same chromosome and have same repeat unit
    var reps: seq[tread] = group.v

    if reps[0].tid < 0:
      #stderr.write_line "yielding " & $reps.len & " unplaced reads with repeat: " & $reps[0].repeat
      yield Cluster(reads: reps)
      continue

    for c in trcluster(reps, max_dist, min_supporting_reads):
      yield c


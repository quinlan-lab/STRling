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

type tid_rep* = tuple[tid:int32, repeat: array[6, char]]

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

type Cluster* = object
  reads*: seq[tread]
  left_most*: uint32
  right_most*: uint32

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
  left_most*:uint32
  right*: uint32
  right_most*: uint32
  center_mass*: uint32
  n_left*: uint16
  n_right*: uint16
  n_total*: uint16
  repeat*: string
  name*: string
  force_report*: bool

const bounds_header* = "#chrom\tleft\tright\trepeat\tname\tleft_most\tright_most\tcenter_mass\tn_left\tn_right\tn_total"

proc `==`(a,b: Bounds): bool =
  if (a.tid == b.tid) and (a.left == b.left) and (a.right == b.right) and (a.repeat == b.repeat):
    return true

# Check if two Bounds overlap. Assumes left <= right in both Bounds
proc overlaps*(a,b: Bounds): bool =
  if a.tid == b.tid and a.repeat == b.repeat:
    var ileft = max(a.left, b.left) #lower bound of intersection interval
    var iright = min(a.right, b.right) #upper bound of intersection interval
    return ileft <= iright #interval non-empty?

# Check if a read (Record) overlaps with a Bounds overlap.
# Same logic/assumptions as overlaps for Bounds above
proc overlaps*(a: Record, b: Bounds): bool =
  if a.tid == b.tid:
    var ileft = max(a.start, int(b.left))
    var iright = min(a.stop, int(b.right))
    return ileft <= iright

# Parse single line of an STR loci bed file
proc parse_bedline*(l:string, targets: seq[Target], window: uint32): Bounds =
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
  if result.repeat.len > 6:
    quit &"ERROR: STRling currently only supports 1-6 bp repeat units. Input bed contains repeat unit length {result.repeat.len}\n{l}"
  # Ensure left_most and right_most are within the chromosome
  # (convert uint to int to allow negative numbers)
  result.left_most = uint32(max(int32(result.left) - int32(window), 0))
  result.right_most = min(result.right + window, targets[result.tid].length)

  for x in result.repeat:
    if x notin "ATCG":
      quit fmt"Error reading loci bed file. Expected DNA (ATCG only) in the 4th field, and got an unexpected character on line: {l}"
  doAssert(result.left <= result.right, &"{result}")
  doAssert(result.left_most <= result.right_most, &"{result}")

# Parse an STR loci bed file
proc parse_bed*(f:string, targets: seq[Target], window: uint32): seq[Bounds] =
  for line in lines f:
    result.add(parse_bedline(string(line), targets, window))

# Parse single line of a STRling bounds file
proc parse_boundsline*(l:string, targets: seq[Target]): Bounds =
  var l_split = l.split("\t")
  if len(l_split) != 11:
    quit fmt"Error reading loci bed file. Expected 11 fields and got {len(l_split)} on line: {l}"
  result.tid = int32(get_tid(l_split[0], targets))
  result.left = uint32(parseInt(l_split[1]))
  result.right = uint32(parseInt(l_split[2]))
  result.repeat = l_split[3]
  result.name = l_split[4]
  result.left_most = uint32(parseInt(l_split[5]))
  result.right_most = uint32(parseInt(l_split[6]))
  result.center_mass = uint32(parseInt(l_split[7]))
  result.n_left = uint16(parseInt(l_split[8]))
  result.n_right = uint16(parseInt(l_split[9]))
  result.n_total = uint16(parseInt(l_split[10]))
  for x in result.repeat:
    if x notin "ATCG":
      quit fmt"Error reading loci bed file. Expected DNA (ATCG only) in the 4th field, and got an unexpected character on line: {l}"
  doAssert(result.left <= result.right, &"{l}")
  doAssert(result.left_most <= result.right_most, &"{l}")

# Parse an STRling bounds file
proc parse_bounds*(f:string, targets: seq[Target]): seq[Bounds] =
  for line in lines f:
    if line.startswith("#"): continue
    result.add(parse_boundsline(string(line), targets))


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

  # Set the positions of the left and right most informative reads
  if int(cl.left_most) > 0:
    result.left_most = cl.left_most
  else:
    result.left_most = posns.min()
  if int(cl.right_most) > 0:
    result.right_most = cl.right_most
  else:
    result.right_most = posns.max()

  # XXX this correction may be hiding a bug elsewhere
  if result.left_most > result.left:
    result.left_most = result.left
  if result.right_most < result.right:
    result.right_most = result.right

  doAssert(result.left <= result.right, &"{result}")
  doAssert(result.left_most <= result.right_most, &"{result}")

proc trim(cl:var Cluster, max_dist:uint32) =
  if cl.reads.len == 0: return
  # drop stuff from start of cluster that is now outside the expected distance
  var lo = max(0, cl.posmed(mediani).int - max_dist.int).uint32
  while len(cl.reads) > 1 and cl.reads[0].position < lo:
    cl.reads = cl.reads[1..cl.reads.high]

proc id*(b:Bounds, targets: seq[Target]): string =
  return &"{targets[b.tid].name}-{b.left}-{b.repeat}"

proc tostring*(b:Bounds, targets: seq[Target]): string =
  doAssert(b.left_most <= b.right_most, &"{b}")
  doAssert(b.left_most <= b.left, &"{b}")
  doAssert(b.right_most >= b.right, &"{b}")
  return &"{targets[b.tid].name}\t{b.left}\t{b.right}\t{b.repeat}\t{b.name}\t{b.left_most}\t{b.right_most}\t{b.center_mass}\t{b.n_left}\t{b.n_right}\t{b.n_total}"

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

      c1.right_most = mid - 1
      c2.left_most = mid

      yield c1
      yield c2

    else:
      yield c


iterator trcluster*(reps: seq[tread], max_dist:uint32, min_supporting_reads:int): Cluster =
  var i = 0
  var c:Cluster
  var last_right = 0'u32
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
      c.right_most = max(c.reads[c.reads.high].position, c.posmed(mediani) + max_dist)
      c.left_most = min(c.reads[0].position, c.posmed(mediani) - max_dist)

      if c.reads.len >= min_supporting_reads and c.reads.has_anchor:
        for sc in c.split_cluster(min_supporting_reads):
          yield sc
          last_right = sc.right_most
        c = Cluster()
      # increment i to past last j and break out of this cluster
      break

  c.trim(max_dist + 100)
  c.right_most = max(c.reads[c.reads.high].position, c.posmed(mediani) + max_dist)
  c.left_most = min(c.reads[0].position, c.posmed(mediani) - max_dist)

  doAssert(c.left_most <= c.right_most)

  if c.reads.len >= min_supporting_reads and c.reads.has_anchor:
    for sc in c.split_cluster(min_supporting_reads):
      yield sc

type tread_id* = object
  tr*: tread
  id*: uint32

proc has_per_sample_reads(c:Cluster, supporting_reads:int, qname2sample:TableRef[string, uint32]): bool =
  ## check that, within the cluster there are at least `supporting reads` from
  ## 1 sample.
  var sample_counts = initCountTable[uint32](8)
  for r in c.reads:
    sample_counts.inc(qname2sample[r.qname])

  return sample_counts.largest().val >= supporting_reads

iterator cluster*(reps: var seq[tread], max_dist:uint32, min_supporting_reads:int=5): Cluster =
  # reps passed here are guaranteed to be split by tid and repeat unit.
  if reps.len > 0:
    doAssert reps[0].tid == reps[reps.high].tid and reps[0].repeat == reps[reps.high].repeat

    if reps[0].tid < 0:
      #stderr.write_line "yielding " & $reps.len & " unplaced reads with repeat: " & $reps[0].repeat
      yield Cluster(reads: reps)
    else:
      for c in trcluster(reps, max_dist, min_supporting_reads):
        yield c

iterator cluster*(id_reps: var seq[tread_id], max_dist:uint32, supporting_reads:int=5): Cluster =
  # reps passed here are guaranteed to be split by tid and repeat unit.
  # need a lookup from qname -> sample and to extract the qreads from each sample.
  var qname2sample = newTable[string, uint32]()
  var reps = newSeqOfCap[tread](id_reps.len)
  for r in id_reps:
    qname2sample[r.tr.qname] = r.id
    reps.add(r.tr)

  # then call normal cluster iterator
  for c in reps.cluster(max_dist, supporting_reads):
    # then check if this cluster has enough reads from 1 sample.
      if c.has_per_sample_reads(supporting_reads, qname2sample):
        yield c

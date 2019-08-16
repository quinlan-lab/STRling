import algorithm
import math
import strformat
import tables
import itertools
import hts/bam
import strutils

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

# Parse a regions file
proc parse_bounds*(l:string): Bounds =
  var l_split = l.splitWhitespace()
  result.tid = int32(parseInt(l_split[0]) - 1) #XXX this will only work for number chomosomes, need to convert this properly
  result.left = uint32(parseInt(l_split[1]))
  result.right = uint32(parseInt(l_split[2]))
  result.repeat = l_split[3]

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
    else:
      posns.add(r.position)
  if posns.len > 0:
    result.center_mass = posns[int(posns.len / 2)]
  if lefts.len > 0:
    var ll = lefts.largest
    result.left = ll.key
  else:
    result.left = result.center_mass
  if rights.len > 0:
    var rr = rights.largest
    result.right = rr.key
  else:
    result.right = result.center_mass

  # If one bound is missing, replace it with the other
  if (result.left == 0) and (result.right > 0'u32):
    result.left = result.right
  if (result.right == 0) and (result.left > 0'u32):
    result.right = result.left

proc trim(cl:var Cluster, max_dist:uint32) =
  if cl.reads.len == 0: return
  # drop stuff from start of cluster that is now outside the expected distance
  var lo = max(0, cl.posmed(mediani).int - max_dist.int).uint32
  while len(cl.reads) > 1 and cl.reads[0].position < lo:
    cl.reads = cl.reads[1..cl.reads.high]

proc tostring*(b:Bounds, targets: seq[Target]): string =
  return &"{targets[b.tid].name}\t{b.left}\t{b.right}\t{b.center_mass}\t{b.n_left}\t{b.n_right}\t{b.n_total}\t{b.repeat}"

proc tostring*(c:Cluster, targets: seq[Target]): string =
  var rep: string
  for v in c.reads[0].repeat:
    if v == 0.char: continue
    rep.add(v)
  return &"{targets[c.reads[0].tid].name}\t{c.reads[0].position}\t{c.reads[c.reads.high].position}\t{c.reads.len}\t{rep}"

iterator cluster*(tandems: var seq[tread], max_dist:uint32, min_supporting_reads:int=5): Cluster =
  tandems.sort(tread_cmp)

  for group in groupby(tandems, bytidrep):
    # reps are on same chromosome and have same repeat unit
    var reps: seq[tread] = group.v

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
        if c.reads.len >= min_supporting_reads:
          yield c
          c = Cluster()
        # increment i to past last j and break out of this cluster
        break

    c.trim(max_dist)
    if c.reads.len >= min_supporting_reads:
      yield c

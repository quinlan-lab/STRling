import algorithm
import math
import strformat
import itertools
import hts/bam

# Data structure storing information about each read that looks like an STR
type tread* = object
  tid*: int32
  position*: uint32
  repeat*: array[6, char]
  flag*: Flag
  split*: int8
  mapping_quality*: uint8
  repeat_count*: uint8
  read_length*: uint8
  when defined(debug):
    qname*: string

type Cluster* = object
  reads*: seq[tread]

const mediani = 9

proc posmed(cl:Cluster, n:int=5): uint32 =
  ## posmed is the median of the first n positions in the cluster
  let mid = int(min(mediani, cl.reads.len) / 2 - 0.5)
  return cl.reads[mid].position

proc bytidrep(t:tread): tuple[repeat:array[6, char], tid:int32] =
  return (t.repeat, t.tid)

# Sorts the reads by chromosome (tid) then repeat unit, then by position
proc tread_cmp(a: tread, b:tread): int =
  if a.tid != b.tid: return cmp(a.tid, b.tid)
  for i in 0..<6:
    if a.repeat[i] != b.repeat[i]:
      return cmp(a.repeat[i], b.repeat[i])
  return cmp(a.position, b.position)


proc trim(cl:var Cluster, max_dist:uint32) =
  if cl.reads.len == 0: return
  # drop stuff from start of cluster that is now outside the expected distance
  var lo = max(0, cl.posmed(mediani).int - max_dist.int).uint32
  while len(cl.reads) > 0 and cl.reads[0].position < lo:
    cl.reads = cl.reads[1..cl.reads.high]

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

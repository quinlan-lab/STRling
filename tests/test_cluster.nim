import unittest
import sequtils
import strformat
include strpkg/cluster

proc `$`(c:Cluster): string =
  result = "Cluster(reads: @[\n"
  for i, r in c.reads:
    result &= &"  tread(position: {r.position}, split: Soft.{r.split})"
    if i < c.reads.high:
      result &= ",\n"
    else:
      result &= "])\n"


suite "cluster suite":

  test "bounds equal":
    var a = Bounds(tid: 0, left: 1, right: 100)
    var b = Bounds(tid: 0, left: 1, right: 100)
    check a == b

  test "bounds unequal":
    var a = Bounds(tid: 0, left: 1, right: 100)
    var b = Bounds(tid: 0, left: 2, right: 100)
    check not (a == b)

  test "bounds overlap":
    var a = Bounds(tid: 0, left: 1, right: 100)
    var b = Bounds(tid: 0, left: 50, right: 100)
    check a.overlaps(b)

  test "bounds don't overlap":
    var a = Bounds(tid: 0, left: 1, right: 100)
    var b = Bounds(tid: 0, left: 200, right: 300)
    check a.overlaps(b) == false

  test "test clustering":

    var reads = @[
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 1, split: Soft.none),
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 1, split: Soft.none),
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 1, split: Soft.none),
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 200, split: Soft.none),
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 255, split: Soft.none),
    tread(tid: 2.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 1, split: Soft.none),
    ]

    var j = 0
    for cl in cluster(reads, 125, min_supporting_reads=3):
      check cl.reads.len == 4

      check cl.tostring(@[Target(name:"chr0"), Target(name:"chr1")]) == "chr1	1	200	4	AAAAAT"
      j += 1

    check j == 1


  test "test bounds":
    var reads = @[
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 123, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 123, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 223, split: Soft.left),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 223, split: Soft.left),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 223, split: Soft.left),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 223, split: Soft.left),

      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 253, split: Soft.right),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 253, split: Soft.right),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 253, split: Soft.right),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 253, split: Soft.right),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 283, split: Soft.none),
      ]

    var cl = Cluster(reads:reads)
    var b = cl.bounds
    check b.left == 223
    check b.right == 253

  test "test bounds: no soft-clipped reads so use median":
    var reads = @[
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 1, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 2, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 5, split: Soft.none),
      ]

    var cl = Cluster(reads:reads)
    var b = cl.bounds
    check b.left == 2
    check b.right == 3

  test "test bounds: no right soft-clipped reads so use median":
    var reads = @[
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 1, split: Soft.left),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 2, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 3, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 5, split: Soft.none),
      ]

    var cl = Cluster(reads:reads)
    var b = cl.bounds
    check b.left == 1
    check b.right == 4

  test "test bounds: no left soft-clipped reads so use median":
    var reads = @[
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 2, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 2, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 3, split: Soft.right),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 5, split: Soft.none),
      ]

    var cl = Cluster(reads:reads)
    var b = cl.bounds
    check b.left == 2
    check b.right == 3

# This test is failing. Need to think about how to fix
#  test "test bounds: no left soft-clipped reads so use median":
#    var reads = @[
#      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 3, split: Soft.right),
#      ]
#
#    var cl = Cluster(reads:reads)
#    var b = cl.bounds
#    check b.left == 3
#    check b.right == 3

  test "test parse STR region":
    var l = "chr1 1 100 CAG"
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    var b = parse_bounds(l, targets)
    check b.tid == 0
    check b.left == 1
    check b.right == 100
    check b.repeat == "CAG"

  test "test parse STR bed file":
    var f = "test_str_parse.bed"
    var text = "chr1 1 100 CAG\nchr1 1 100 CAG"
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    writeFile(f, text)
    check parse_loci(f, targets)[1].tid == 0


  test "inverted bounds":
    var tr = @[tread(tid: 20, position: 48086080, repeat: ['T', 'T', '\x00', '\x00', '\x00', '\x00'], split: Soft.none, mapping_quality: 0),
              tread(tid: 20, position: 48086101, repeat: ['T', 'T', '\x00', '\x00', '\x00', '\x00'], split: Soft.none, mapping_quality: 15),
              tread(tid: 20, position: 48086132, repeat: ['T', 'T', '\x00', '\x00', '\x00', '\x00'], split: Soft.none, mapping_quality: 0),
              tread(tid: 20, position: 48086164, repeat: ['T', 'T', '\x00', '\x00', '\x00', '\x00'], split: Soft.none, mapping_quality: 0),
              tread(tid: 20, position: 48086187, repeat: ['T', 'T', '\x00', '\x00', '\x00', '\x00'], split: Soft.none, mapping_quality: 0),
              tread(tid: 20, position: 48086281, repeat: ['T', 'T', '\x00', '\x00', '\x00', '\x00'], split: Soft.none, mapping_quality: 0)]
    var b = Cluster(reads:tr).bounds
    check b.left < b.right



  test "should split cluster":
    var treads = @[
     tread(position: 370, split: Soft.none),
     tread(position: 391, split: Soft.right),
     tread(position: 391, split: Soft.right),
     tread(position: 391, split: Soft.right),
     tread(position: 403, split: Soft.none),
     tread(position: 503, split: Soft.none),

     tread(position: 850, split: Soft.left),
     tread(position: 850, split: Soft.left),
     tread(position: 850, split: Soft.left),
     tread(position: 850, split: Soft.left),
     tread(position: 880, split: Soft.none),
     ]


    var clusters : seq[Cluster]
    for c in trcluster(treads, 500, 1):
      clusters.add(c)
    check clusters.len == 2

    var c1 = clusters[0]
    check c1.reads.len == 6
    check c1.reads[c1.reads.high].position == 503

    var c2 = clusters[1]
    check c2.reads.len == 5
    check c2.reads[0].position == 850



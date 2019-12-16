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
    check b.left_most == 123
    check b.right_most == 283

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
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 1, split: Soft.left),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 2, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 3, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 5, split: Soft.none),
      ]

    var cl = Cluster(reads:reads)
    var b = cl.bounds
    check b.left == 1
    check b.right == 2

  test "test bounds: no left soft-clipped reads so use median":
    var reads = @[
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 2, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 2, split: Soft.none),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 3, split: Soft.right),
      tread(tid:1'i32, repeat: ['A', 'T', 'G', '\x00', '\x00', '\x00'], position: 5, split: Soft.none),
      ]

    var cl = Cluster(reads:reads)
    var b = cl.bounds
    check b.left == 3
    check b.right == 4

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

  test "test parse line from STR bed file":
    var l = "chr1 100 200 CAG"
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    var window = 50
    var b = parse_bedline(l, targets, window.uint32)
    check b.tid == 0
    check b.left == 100
    check b.left_most == 50
    check b.right == 200
    check b.right_most == 250
    check b.repeat == "CAG"

  test "test parse STR bed file":
    var f = "test_str_parse.bed"
    var text = "chr1 1 100 CAG\nchr1 1 100 CAG"
    var window = 100
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    writeFile(f, text)
    check parse_bed(f, targets, window.uint32)[1].tid == 0

  test "test parse line from STRling -bounds.txt":
    var l = "chr1\t990\t1010\tCAG\t\t500\t1500\t1000\t3\t1\t50"
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    var b = parse_boundsline(l, targets)
    check b.tid == 0
    check b.left == 990
    check b.right == 1010
    check b.repeat == "CAG"

  test "test parse STRling -bounds.txt file":
    var f = "test_str_parse-bounds.txt"
    var text = "chr1\t990\t1010\tCAG\t\t500\t1500\t1000\t3\t1\t50\nchr1\t990\t1010\tCAG\t\t500\t1500\t1000\t3\t1\t50"
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    writeFile(f, text)
    check parse_bounds(f, targets)[1].tid == 0


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


  test "inverted bounds again":

    var tr = @[
      tread(tid: 11, position: 115977335, split: Soft.none),
      tread(tid: 11, position: 115977397, split: Soft.none),
      tread(tid: 11, position: 115977419, split: Soft.none),
      tread(tid: 11, position: 115977448, split: Soft.left),
      tread(tid: 11, position: 115977585, split: Soft.none),
      tread(tid: 11, position: 115977598, split: Soft.none)]

    var b = Cluster(reads: tr).bounds
    check b.left < b.right

  test "inverted bounds 3":
    var tr = @[
       tread(tid: 10, position: 92611809, split: Soft.none),
       tread(tid: 10, position: 92611833, split: Soft.right),
       tread(tid: 10, position: 92611833, split: Soft.right),
       tread(tid: 10, position: 92611921, split: Soft.none),
       tread(tid: 10, position: 92611939, split: Soft.none)]
    var b = Cluster(reads: tr).bounds
    check b.left < b.right

  test "right_most bug":
    var tr = @[
      tread(tid: 5, position: 34847227, split: Soft.left),
      tread(tid: 5, position: 34847227, split: Soft.none),
      tread(tid: 5, position: 34847883, split: Soft.left),
      tread(tid: 5, position: 34847911, split: Soft.none),
      tread(tid: 5, position: 34847921, split: Soft.left),
      tread(tid: 5, position: 34847921, split: Soft.left),
      tread(tid: 5, position: 34847930, split: Soft.none),
      tread(tid: 5, position: 34848950, split: Soft.left),
      tread(tid: 5, position: 34848950, split: Soft.left),
      tread(tid: 5, position: 34848950, split: Soft.left)]

    var b = Cluster(reads: tr).bounds
    check b.left < b.right

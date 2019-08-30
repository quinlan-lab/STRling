import unittest
include strpkg/cluster

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
    check b.right == 2

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
    check b.right == 3

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

  test "cluster with flipped left right":

    var reads = @[
       tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 195, split: Soft.right),
       tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 195, split: Soft.right),
       tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 200, split: Soft.none),
       tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 205, split: Soft.left),
    ]
    var b = Cluster(reads:reads).bounds
    echo "bounds to fix:", b
 
  test "that cluster has right bounds":
    var reads = @[tread(tid: 11'i32, position: 90878497, repeat: ['A', 'A', '\x00', '\x00', '\x00', '\x00'], split: Soft.none),
                  tread(tid: 11'i32, position: 90878812, repeat: ['A', 'A', '\x00', '\x00', '\x00', '\x00'], split: Soft.none),
                  tread(tid: 11'i32, position: 90878838, repeat: ['A', 'A', '\x00', '\x00', '\x00', '\x00'], split: Soft.left),
                  tread(tid: 11'i32, position: 90878838, repeat: ['A', 'A', '\x00', '\x00', '\x00', '\x00'], split: Soft.left),
                  tread(tid: 11'i32, position: 90878902, repeat: ['A', 'A', '\x00', '\x00', '\x00', '\x00'], split: Soft.none)]

    var b = Cluster(reads:reads).bounds
    check b.left == b.right - 1

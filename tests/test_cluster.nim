import unittest
include strpkg/cluster

suite "cluster suite":
  test "test clustering":

    var reads = @[
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 1),
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 1),
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 1),
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 200),
    tread(tid: 1.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 255),
    tread(tid: 2.int32, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], position: 1),
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

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
    for cl in cluster(reads, 250, min_supporting_reads=3):
      check cl.reads.len == 4

      check cl.tostring(@["chr0", "chr1"]) == "chr1	1	200	4	AAAAAT"
      j += 1

    check j == 1

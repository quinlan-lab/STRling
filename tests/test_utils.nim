import unittest
include strpkg/utils

suite "utils suite":

  test "test get_tid":
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    check get_tid("chr1", targets) == 0



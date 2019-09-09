import unittest
include strpkg/utils

suite "utils suite":

  test "test get_tid":
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    check get_tid("chr1", targets) == 0

  test "test median":
    check median_depth(@[1,2,2]) == 2
    check median_depth(@[2000]) == 1047 #Values greater than 1047 will be set to 1047 to avoid overflow
    #echo median_depth(@[-11,2,2]) # This will throw an inded out of range error

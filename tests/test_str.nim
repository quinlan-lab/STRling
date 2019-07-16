import unittest
include str
import hts/bam

suite "str suite":
  test "debug has qname":
    when defined(debug):
      var tr = tread(qname:"my_qname", tid:0'i32, position: 222'u32, split: 0'i8, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], repeat_count: 123'u8)
      check tr.tostring(@[Target(name:"chr1")]) == "chr1	222	AAAAAT	0	123	my_qname"


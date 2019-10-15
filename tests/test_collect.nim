import unittest
import strpkg/collect
import strpkg/cluster
import hts/bam

suite "collect suite":

  test "test spanning read":
    
    var h = Header()
    h.from_string("""@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000""")
    check h.hdr != nil

    var a = NewRecord(h)
    var txt = "read1	0	chr1	1	40	15M5S	*	0	0	AAAAAAAAAAAAAAAAAAAA	*"
    a.from_string(txt)

    var support:Support
    var bounds:Bounds

    bounds = Bounds(tid: 0, left: 5, right: 100, repeat: "A")
    check spanning_read(a, bounds, support) == false

    bounds = Bounds(tid: 0, left: 5, right: 10, repeat: "A")
    check spanning_read(a, bounds, support) == true

    bounds = Bounds(tid: 0, left: 1, right: 15, repeat: "A")
    check spanning_read(a, bounds, support) == false

    # Require the read to span (repeat length - 1) bp either side of the bounds
    bounds = Bounds(tid: 0, left: 5, right: 7, repeat: "AAAAAA")
    check spanning_read(a, bounds, support) == false

    bounds = Bounds(tid: 0, left: 6, right: 7, repeat: "AAAAAA")
    check spanning_read(a, bounds, support) == true


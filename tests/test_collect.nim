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

  test "test spanning pair":
    
    var h = Header()
    h.from_string("""@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:10000""")
    check h.hdr != nil

    var a = NewRecord(h)
    a.from_string("read1	99	chr1	1	40	15M5S	chr1	500	0	AAAAAAAAAAAAAAAAAAAA	*")
    var b = NewRecord(h)
    b.from_string("read1	147	chr1	500	40	15M5S	chr1	1	0	AAAAAAAAAAAAAAAAAAAA	*")
    
    var frag_sizes: array[4096, uint32]
    frag_sizes[0] = 500
    frag_sizes[1] = 450
    var support:Support
    var bounds:Bounds
    bounds = Bounds(tid: 0, left: 100, right: 150, repeat: "A")

    check spanning_fragment(a, b, bounds, support, frag_sizes) == true



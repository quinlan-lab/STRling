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
    var txt = "read1	0	chr1	1	40	25M5S	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*"
    a.from_string(txt)

    var support:Support
    var bounds:Bounds

    bounds = Bounds(tid: 0, left: 5, right: 100, repeat: "A")
    check spanning_read(a, bounds, support) == false

    # Require the read to span (repeat length - 1) bp either side of the bounds
    bounds = Bounds(tid: 0, left: 5, right: 15, repeat: "AAAAAA")
    check spanning_read(a, bounds, support) == false

    bounds = Bounds(tid: 0, left: 6, right: 15, repeat: "AAAAAA")
    check spanning_read(a, bounds, support) == true

    # Require extra slop for small bounds
    bounds = Bounds(tid: 0, left: 9, right: 10, repeat: "AAAAAA")
    check spanning_read(a, bounds, support) == false

    bounds = Bounds(tid: 0, left: 10, right: 11, repeat: "AAAAAA")
    check spanning_read(a, bounds, support) == true

  test "test spanning pair":
    
    var h = Header()
    h.from_string("""@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:10000""")
    check h.hdr != nil

    var a = NewRecord(h)
    var atxt = "read1	99	chr1	1	40	15M5S	chr1	500	0	AAAAAAAAAAAAAAAAAAAA	*"
    a.from_string(atxt)
    var b = NewRecord(h)
    var btxt = "read1	147	chr1	500	40	15M5S	chr1	1	0	AAAAAAAAAAAAAAAAAAAA	*"
    b.from_string(btxt)
    
    var frag_sizes: array[4096, uint32]
    var support:Support
    var bounds:Bounds

    bounds = Bounds(tid: 0, left: 100, right: 150, repeat: "A")
    check spanning_fragment(a, b, bounds, support, frag_sizes) == true

    bounds = Bounds(tid: 0, left: 450, right: 513, repeat: "A")
    check spanning_fragment(a, b, bounds, support, frag_sizes) == true

    # Require extra slop for small bounds
    bounds = Bounds(tid: 0, left: 512, right: 513, repeat: "A")
    check spanning_fragment(a, b, bounds, support, frag_sizes) == false



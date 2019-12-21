import unittest
import strpkg/collect
import strpkg/cluster
import hts/bam

suite "collect suite":

  test "test overlapping and spanning reads":
    
    var h = Header()
    h.from_string("""@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000""")
    check h.hdr != nil

    var a = NewRecord(h)
    var txt = "read1	0	chr1	1	40	25M5S	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*"
    a.from_string(txt)

    var bounds:Bounds

    var s1:Support
    bounds = Bounds(tid: 0, left: 50, right: 100, repeat: "A")
    check overlapping_read(a, bounds, s1) == false

    # Require the read to span (repeat length - 1) bp either side of the bounds
    var s2:Support
    bounds = Bounds(tid: 0, left: 5, right: 15, repeat: "AAAAAA")
    check overlapping_read(a, bounds, s2) == true
    check s2.Type == SupportType.OverlappingRead

    var s3:Support
    bounds = Bounds(tid: 0, left: 6, right: 15, repeat: "AAAAAA")
    check overlapping_read(a, bounds, s3) == true
    check s3.Type == SupportType.SpanningRead

    # Require extra slop for small bounds
    var s4:Support
    bounds = Bounds(tid: 0, left: 9, right: 10, repeat: "AAAAAA")
    check overlapping_read(a, bounds, s4) == true
    check s2.Type == SupportType.OverlappingRead

    var s5:Support
    bounds = Bounds(tid: 0, left: 10, right: 11, repeat: "AAAAAA")
    check overlapping_read(a, bounds, s5) == true
    check s5.Type == SupportType.SpanningRead

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



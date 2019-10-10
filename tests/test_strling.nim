import unittest
import tables
import strling
import strpkg/utils
import strpkg/extract
import strpkg/call
import hts/bam

suite "strling suite":
  test "debug has qname":
    when defined(debug):
      var tr = tread(qname:"my_qname", tid:0'i32, position: 222'u32, split: Soft.none, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], repeat_count: 123'u8)
      check tr.tostring(@[Target(name:"chr1")]) == "chr1	222	AAAAAT	none	123	my_qname"


  test "reset count":
    var h = Header()

    h.from_string("""@HD	VN:1.3	SO:coordinate
@SQ	SN:X	LN:249250621
@SQ	SN:18	LN:249250621""")
    check h.hdr != nil
    var a = NewRecord(h)

    var txt = "1-2014	97	18	5972849	42	20S127M3S	X	9760012	0	AATAGAATAGAATAAAATAGAATAGAATAGTATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAACAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAG	=CCGGGGGGGGGG=(GJJJJJ=JJCJJJJJ1JCGGJGJJJGGGJJGJGJGGJGJGJGGCCGCJGGCGGJJGCGGGGCGG=CGGGGGG(GGCCGC=CCGGCJCCGGGGGCCGGGGCGGGGGCGG1G8CGG=CGGGG==GGCGCGGGG=GGC	NM:i:3	MD:Z:37C4G24T59	MC:Z:150M	AS:i:112	XS:i:92	RG:Z:disease_loci_sims_minpath.sorted.22-46191235_ATTCT_0_850.22-46191235_0.5_1.stutter.merged_L001	XA:Z:2,-61912014,40S57M1I52M,3;"
    a.from_string(txt)
    echo a

    var cache = Cache(tbl:newTable[string, tread](8192), cache: newSeqOfCap[tread](65556))
    var opts = Options(median_fragment_length: 500, proportion_repeat: 0.6, min_mapq: 20'u8)
    var counts = init[uint8]()
    var rep: array[6, char]
    rep[0] = 'A'

    cache.add_soft(a, counts, opts, rep)
    check cache.cache.len == 1

    txt = "1-2014	145	X	9760012	60	150M	18	5972849	0	AATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAGAATAG	CGCCGCGG=GGGCGGC=GCCCCGGGGGGGCC=G=CC===CGCGJCJ=J=8GGCGGGGGCGGCGCGC8GGGCGGGGGC8CGGJGGGGG=JJJGCGG8JGJ1JJJGGJGJGGJGGGGJJJCJJJCGJ8JJGJJJJJJJJCG=CGGGGGG1CC	NM:i:0	MD:Z:150	MC:Z:20S127M3S	AS:i:150	XS:i:32	RG:Z:disease_loci_sims_minpath.sorted.22-46191235_ATTCT_0_850.22-46191235_0.5_1.stutter.merged_L001"

    a.from_string(txt)
    cache.add_soft(a, counts, opts, rep)
    check cache.cache.len == 1


  test "monomer repeat":
    var h = Header()

    h.from_string("""@HD	VN:1.3	SO:coordinate
@SQ	SN:X	LN:249250621
@SQ	SN:18	LN:249250621""")
    check h.hdr != nil
    var a = NewRecord(h)

    var txt = "1-2014	97	18	5972849	42	20S127M3S	X	9760012	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	=CCGGGGGGGGGG=(GJJJJJ=JJCJJJJJ1JCGGJGJJJGGGJJGJGJGGJGJGJGGCCGCJGGCGGJJGCGGGGCGG=CGGGGGG(GGCCGC=CCGGCJCCGGGGGCCGGGGCGGGGGCGG1G8CGG=CGGGG==GGCGCGGGG=GGC	NM:i:3	MD:Z:37C4G24T59	MC:Z:150M	AS:i:112	XS:i:92	RG:Z:disease_loci_sims_minpath.sorted.22-46191235_ATTCT_0_850.22-46191235_0.5_1.stutter.merged_L001	XA:Z:2,-61912014,40S57M1I52M,3;"
    a.from_string(txt)

    var cache = Cache(tbl:newTable[string, tread](8192), cache: newSeqOfCap[tread](65556))
    var opts = Options(median_fragment_length: 500, proportion_repeat: 0.6, min_mapq: 20'u8)
    var counts = init[uint8]()
    var repeat_count = 0
    var align_length = 0

    var rep = a.get_repeat(nil, counts, repeat_count, align_length, opts)
    check rep == ['A', 'A', '\x00', '\x00', '\x00', '\x00']


  test "unplaced read pair: both str":
    var opts = Options(median_fragment_length: 500, proportion_repeat: 0.8, min_mapq: 20'u8)
    var A = tread(qname:"my_qname", tid:0'i32, position: 222'u32, split: Soft.none, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], repeat_count: 150'u8, mapping_quality: 30)
    var B = tread(qname:"my_qname", tid:0'i32, position: 222'u32, split: Soft.none, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], repeat_count: 150'u8, mapping_quality: 30)
    check unplaced_pair(A, B, opts) == true 

  test "unplaced read pair: one str, one low mapq":
    var opts = Options(median_fragment_length: 500, proportion_repeat: 0.8, min_mapq: 20'u8)
    var A = tread(qname:"my_qname", tid:0'i32, position: 222'u32, split: Soft.none, repeat: ['A', 'A', 'A', 'A', 'A', 'T'], repeat_count: 150'u8, mapping_quality: 16)
    var B = tread(qname:"my_qname", tid:0'i32, position: 222'u32, split: Soft.none, repeat: ['\x00', '\x00', '\x00', '\x00', '\x00', '\x00'], repeat_count: 0'u8, mapping_quality: 16)
    check unplaced_pair(A, B, opts) == true 

  test "unplaced read pair: false":
    var opts = Options(median_fragment_length: 500, proportion_repeat: 0.8, min_mapq: 20'u8)
    var A = tread(qname:"my_qname", tid:0'i32, position: 222'u32, split: Soft.none, repeat: ['\x00', '\x00', '\x00', '\x00', '\x00', '\x00'], repeat_count: 150'u8, mapping_quality: 30)
    var B = tread(qname:"my_qname", tid:0'i32, position: 222'u32, split: Soft.none, repeat: ['\x00', '\x00', '\x00', '\x00', '\x00', '\x00'], repeat_count: 0'u8, mapping_quality: 30)
    check unplaced_pair(A, B, opts) == false


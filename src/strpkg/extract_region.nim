import hts/bam
import argparse
import algorithm
import strformat
import tables

proc get_mate(aln:Record, ibam:Bam): Record =
  if aln.mate_tid == -1:
    for o in ibam.query("*"):
      if o.flag.supplementary or o.flag.secondary: continue
      if o.flag.read1 == aln.flag.read1: continue
      if o.qname == aln.qname:
        return o.copy()
  else:
    for o in ibam.query(aln.mate_tid, aln.mate_pos.int - 1, aln.mate_pos.int + 1):
      if o.flag.supplementary or o.flag.secondary: continue
      if o.flag.read1 == aln.flag.read1: continue
      if o.qname == aln.qname:
        return o.copy()
  stderr.write_line "skipping pair. mate not found for " & aln.tostring

proc extract_region_main*() =

  var p = newParser("strling pull_region"):
    option("-f", "--fasta", help="path to fasta file, only required for cram")
    option("-o", "--output-bam", help="path to output bam", default="extracted.bam")
    arg("bam")
    arg("region")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "pull_region":
    argv = argv[1..argv.high]
  if len(argv) == 0: argv = @["-h"]

  var opts = p.parse(argv)
  if opts.help:
    quit 0

  var ibam:Bam
  if not ibam.open(cstring(opts.bam), fai=cstring(opts.fasta), threads=3, index=true):
    quit "could not open bam"

  var records = newSeq[Record]()
  var counts = newTable[string, int]()

  for aln in ibam.query(opts.region):
    if aln.flag.secondary or aln.flag.supplementary: continue
    records.add(aln.copy())

    counts.mgetOrPut(aln.qname, 0).inc
  stderr.write_line &"extracted {records.len} alignments. now checking for mates"

  var mates: seq[Record]
  for i, r in records:
    if i mod 10000 == 0:
      stderr.write_line &"extracting mates. on records {i} of {records.len}"
    if counts[r.qname] == 2: continue # got both mates

    var m = r.get_mate(ibam)
    if m != nil:
      mates.add(m)

  records.add(mates)

  records.sort(proc(a, b:Record): int =
    if a.tid != b.tid: return cmp(a.tid, b.tid)
    return cmp(a.start, b.start)
  )

  var obam:Bam

  if not obam.open(cstring(opts.output_bam), mode="wb", threads=2):
    quit "couldn't open output bam"

  obam.write_header(ibam.hdr)

  for r in records:
    obam.write(r)

  obam.close

when isMainModule:
  extract_region_main()

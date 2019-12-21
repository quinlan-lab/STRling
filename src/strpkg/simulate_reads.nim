import argparse
import osproc
import hts/fai
import hts/bam
import strutils
import strformat
import random
import ./utils

type Alleles = object
  chrom: string
  position: int
  counts: array[2, int]
  repeat_unit: string

proc parse_allele(str: string): Alleles =
  var toks = str.split(':')
  doAssert toks.len == 3, "error incorrect allele format:" & str
  result.chrom = toks[0]
  result.position = parseInt(toks[1])
  toks = toks[2].split('_')
  doAssert toks.len == 2, "error incorrect allele format:" & str
  result.repeat_unit = toks[0]
  toks = toks[1].split('/')
  result.counts[0] = parseInt(toks[0])
  result.counts[1] = parseInt(toks[1])

randomize()

proc simulate_allele(fai:Fai, fh1:File, fh2:File, allele:var Alleles, fragment_dist: array[4096, uint32], flank: int, depth:int, read_length:int=150) =
  var reference = fai.get(allele.chrom, max(0, allele.position - flank), allele.position + flank + fragment_dist.len)
  var off = reference.find(allele.repeat_unit, start=flank-1, last=flank+2*(1+allele.repeat_unit.len))
  if off == -1:
    off = reference.find(allele.repeat_unit.reverse_complement, start=flank-1, last=flank+2*(1+allele.repeat_unit.len))
    if off == -1:
      stderr.write_line &"warning: couldn't find {allele.repeat_unit} around {allele.chrom}:{allele.position}"
      stderr.write_line &"complement: {allele.repeat_unit.reverse_complement}"
      off = flank
    else:
      allele.repeat_unit = allele.repeat_unit.reverse_complement

  var haplotypes = [reference, reference]
  for i, c in allele.counts:
    var H = haplotypes[i]
    if c == 0: continue
    if c > 0:
      var rep = allele.repeat_unit.repeat(c)
      haplotypes[i] = H[0..<off] & rep & H[off..H.high]
    elif c < 0:
      var rep = allele.repeat_unit.repeat(abs(c))
      if reference.find(rep, start=off) == -1:
        stderr.write_line &"couln't find {c} units of {allele.repeat_unit} around {allele.chrom}:{allele.position} to remove"
      else:
        haplotypes[i] = H[0..<off] & H[(off + rep.len)..H.high]

  # subtract 2 * fragment_dist.len because we added it to avoid going off end.
  #let L = reference.len - 2 * fragment_dist.len
  let L = max(haplotypes[0].len, haplotypes[1].len) - 2 * fragment_dist.len
  let n_total_reads = int(depth.float64 * L.float64 / read_length.float64)
  let n_fragments = int(n_total_reads / 2)

  # convert the fragment distribution back into the raw counts (so in
  # fragment_dist, a value of 233 at index 400 means 233 reads had a fragment
  # length of 400. conver that into 233 entries of 400
  var raw_counts = newSeqOfCap[uint16](65536)
  for i, n in fragment_dist:
    for j in 0..<n:
      raw_counts.add(i.uint16)

  let
     left = 0 #allele.position - flank
     right = max(haplotypes[0].len, haplotypes[1].len) - 2 * fragment_dist.len


  var ihaps = [0, 1]
  var qual = "I".repeat(read_length)

  for i in 0..<n_fragments:
    let fragment_length = sample(raw_counts).int
    var r2_start = -1
    var r1_start = -1
    while r2_start < 0:
       r1_start = rand(left..right).int
       r2_start = r1_start.int + fragment_length - read_length

    let ihap = sample(ihaps)
    let hap = haplotypes[ihap]

    let r1 = hap[r1_start..<r1_start + read_length]
    let r2 = hap[r2_start..<r2_start + read_length].reverse_complement

    let qname = &"@{r1_start + allele.position}_{r2_start + allele.position}_{i}_{ihap}"
    # just re-use ihaps to flip a coin here.
    if sample(ihaps) == 0:
      fh1.write(&"{qname}\n{r1}\n+\n{qual}\n")
      fh2.write(&"{qname}\n{r2}\n+\n{qual}\n")
    else:
      fh2.write(&"{qname}\n{r1}\n+\n{qual}\n")
      fh1.write(&"{qname}\n{r2}\n+\n{qual}\n")

proc read_hist(path: string): array[4096, uint32] =
  #fragment_dist = args.bam_or_hist.read_hist
  var i = 0
  for l in path.lines:
    if i > result.high:
      return
    result[i] = l.strip.parseInt.uint32
    i.inc

proc write_hist(H: array[4096, uint32], path:string) =
  var fh:File
  if not open(fh, path, fmWrite):
    quit "couldn't open output file:" & path
  for v in H:
    fh.write_line($v)
  fh.close

proc simulate_main() =

  var p = newParser("simulate reads"):
    option("--fasta", help="path to fasta file indexed with bwa mem")
    option("--flank", help="distance around event to simulate", default="20000")
    option("--depth", help="integer depth to target", default="30")
    option("--output", help="prefix for output files")
    arg("bam_or_hist", help="bam (or .hist) for fragment-length distribution")
    arg("allele", nargs= -1, help="format: {chrom}:{start}:{repeat_unit}_{count1}/{count2} where counts are losses (-) or gains relative to reference e.g.: 4:41747993:CTG_5/11\nor if this ends in .bed, it's a file with format chrom\tstart\tstop\tCTG_5/11")

  let args = p.parse()

  var
    fai:Fai
    ibam:Bam

  if not fai.open(args.fasta):
    quit "couldn't open fasta"

  var fragment_dist: array[4096, uint32]

  if args.bam_or_hist.endsWith(".hist"):
    fragment_dist = args.bam_or_hist.read_hist
  else:

    if not ibam.open(args.bam_or_hist, fai=args.fasta):
      quit "couldn't open bam"
    var cram_opts = 8191 - SAM_RNAME.int - SAM_RGAUX.int - SAM_QUAL.int - SAM_SEQ.int
    discard ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, cram_opts)
    fragment_dist = ibam.fragment_length_distribution()
    ibam.close
    fragment_dist.write_hist(args.output & ".hist")

  var alleles : seq[Alleles]

  for a in args.allele:
    if a.endsWith(".bed"):
      for l in a.lines:
        if l[0] == '#': continue
        let toks = l.strip.split('\t')
        alleles.add((&"{toks[0]}:{toks[1]}:{toks[3]}").parse_allele)
    else:
      alleles.add(a.parse_allele)

  let flank = parseInt(args.flank)
  let depth = parseInt(args.depth)

  var
    fh1:File
    fh2:File
  if not fh1.open(args.output & "_r1.fastq", fmWrite):
    quit "couldn't open output file"
  if not fh2.open(args.output & "_r2.fastq", fmWrite):
    quit "couldn't open output file"

  for allele in alleles.mitems:
    fai.simulate_allele(fh1, fh2, allele, fragment_dist, flank, depth)

  fh1.close
  fh2.close
  let cmd = &"""bwa mem -R "@RG\tID:sim\tSM:sim" {args.fasta} {args.output}_r1.fastq {args.output}_r2.fastq | samtools sort -o {args.output}.bam && samtools index {args.output}.bam"""
  stderr.write_line execCmd(&"""bash -ec 'set -eo pipefail; {cmd}'""")

when isMainModule:
  simulate_main()




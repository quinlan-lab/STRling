import hts/fai
import random
import argparse
import tables
import lapper
import algorithm
import strutils
import ./cluster
import ./utils
import strformat
import ./read_bed
export read_bed
import hts/fai

type Window = object
  chrom: string
  start: int
  stop: int
  repeat: string


proc trim(w:var Window, dna:string): Window =
  doAssert dna.len == w.stop - w.start
  var ow = w

  var expected:uint64
  for v in w.repeat.slide_by(w.repeat.len):
    expected = v
    break

  # NOTE: this stops trimming on first mis-matching kmer. could require 2 mismatches in a row
  # or something, but this actually works well AFAICT

  # trim left
  for enc in dna.slide_by(w.repeat.len):
    if expected != enc: w.start += w.repeat.len
    else: break

  doAssert w.start < w.stop, &"repeat {w.repeat} not found in expected region for {ow}, {dna}, {w}"

  # trim right
  var dnar = newString(dna.len)
  # reverse the repeat and the dna sequence.
  for i, x in dna:
    # high == len - 1
    dnar[dna.high-i] = x

  var rep = w.repeat
  rep.reverse
  for v in rep.slide_by(w.repeat.len):
    expected = v
    break

  for enc in dnar.slide_by(w.repeat.len):
    if expected != enc: w.stop -= w.repeat.len
    else: break

  doAssert w.start < w.stop, &"repeat {w.repeat} not found in expected region for {ow}, {dna}"
  return w

iterator repeat_windows(fai:Fai, window_size:int, step:int, opts:Options): Window =
  var counts = init[uint8]()
  var repeat_count:int
  for tid in 0..<fai.len:
    var chrom = fai[tid]
    var start = 0
    var L = fai.chrom_len(chrom)
    if L > 2_000_000:
      stderr.write_line &"[strling] finding STR regions on reference chromosome: {chrom}"
    var last_w = Window(stop: -1)
    var chrom_seq = fai.get(chrom).toUpperAscii
    shallow(chrom_seq)
    while start < L:
      var dna = chrom_seq[start..<min(L, start + window_size)]
      var rep = dna.get_repeat(counts, repeat_count, opts)
      if repeat_count > 0:
        var w = Window(chrom: chrom, start: start, stop: start + dna.len, repeat: rep.tostring)
        # allow skipping 1 window.
        if last_w.repeat != w.repeat or w.start > last_w.stop + (window_size-step):
          if last_w.stop != -1 and last_w.stop - last_w.start >= (window_size - step):
            last_w.start = max(0, last_w.start - window_size)
            last_w.stop = min(last_w.stop + window_size, chrom_seq.len)
            yield last_w.trim(chrom_seq[last_w.start..<last_w.stop])
          last_w = w
        else:
          last_w.stop = w.stop

      start += step

    if last_w.stop != -1 and last_w.stop - last_w.start >= (window_size - step):
      last_w.start = max(0, last_w.start - window_size)
      last_w.stop = min(last_w.stop + window_size, chrom_seq.len)
      yield last_w.trim(chrom_seq[last_w.start..<last_w.stop])

proc getTempFile*(): string =
  ## from Nimble
  var tmpdir: string
  if existsEnv("TMPDIR") and existsEnv("USER"):
    tmpdir = joinPath(getEnv("TMPDIR"), getEnv("USER"))
  else:
    tmpdir = getTempDir()
  discard existsOrCreateDir(tmpdir)
  randomize()

  result = tmpdir / $rand(int.high) & $rand(int.high) & ".bed"


proc genome_repeats*(fai:Fai, opts:Options, bed_path:string): TableRef[string, Lapper[region]] =
  ## find repeats in the genome. this helps as we can for example skip
  ## reads that are a perfect match to the genome when that part of the genome
  ## is not a repeat.
  var bed_path = bed_path
  var isTmp = bed_path == ""
  if isTmp:
    bed_path = getTempFile()


  if not fileExists(bed_path):
    var fh:File
    if not fh.open(bed_path, fmWrite):
      quit &"[strling] couldn't open bed file: {bed_path} for writing"
    let window_size: int = 100
    let step: int = 60
    var n = 0
    for w in fai.repeat_windows(window_size, step, opts):
      fh.write_line(&"{w.chrom}\t{w.start}\t{w.stop}\t{w.repeat}")
      n += 1
    fh.close
    stderr.write_line &"[strling] found {n} STR-like regions in the genome"
  else:
    stderr.write_line &"[strling] using existing file {bed_path} for genome repeats"

  result = read_bed(bed_path)
  stderr.write_line &"[strling] got STR repeats from genome into an interval tree"
  if isTmp:
    removeFile(bed_path)

proc genome_repeats*(fasta:string, opts:Options, bed_path:string): TableRef[string, Lapper[region]] =
  var fai:Fai
  if not fai.open(fasta):
    quit &"[strling] couldn't open fasta {fasta} make sure file is present and has a .fai index"
  result = fai.genome_repeats(opts, bed_path)


proc check_reference_size*(b:Bounds, fai:Fai, chrom: string): int =
  var reference = fai.get(chrom, b.left.int, b.right.int + 200)
  var right = int(b.right - b.left)
  let  span = right
  var left = 0
  var nmiss = 0
  var ntotal_miss = 0
  while ntotal_miss < 2 and right < reference.len:
    # allow 2 total misses, but not more than 1 in a row.
    if reference[left..right].count(b.repeat).float < 0.5 * span.float:
      nmiss += 1
      ntotal_miss += 1
      if nmiss == 2:
        right -= span
        break
      left += span
      right += span
      continue
    nmiss = 0
    left += span
    right += span
  right -= span

  result = reference[0..<min(right, reference.len)].count(b.repeat)

proc genome_main*() =

  var p = newParser("str genome-sites"):
    option("-p", "--proportion-repeat", help="proportion of read that is repetitive to be considered as STR", default="0.65")
    arg("fasta", help="path to fasta file")
    arg("bed", help="path to output bed file to be created")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "genome-sites":
    argv = argv[1..argv.high]
  if len(argv) == 0: argv = @["--help"]

  var args = p.parse(argv)
  if args.help:
    quit 0

  var fai:Fai
  if not fai.open(args.fasta):
    quit &"[strling] couldn't open fasta {args.fasta} make sure file is present and has a .fai index"

  var opts = Options(proportion_repeat: parseFloat(args.proportion_repeat))
  discard fai.genome_repeats(opts, args.bed)

when isMainModule:

  #[
  import unittest
  suite "genome strs":
    test "bug":
      var w = Window(chrom: "MT", start: 16460, stop: 16569, repeat: "CACGAT")
      var dna = "ATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGTCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATG"
      echo w.trim(dna)
  ]#

  genome_main()



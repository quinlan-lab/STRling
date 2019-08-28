import cluster
import collect
import utils
import tables
import math
import strformat
import hts/bam

type Event = enum
  Insertion
  Deletion

type Evidence = object
  class: string # class of support e.g. "spanning reads"
  repeat: string
  # allele sizes in bp relative to the reference
  allele1_bp: float
  allele2_bp: float
  # total allele sizes in repeat units
  allele1_ru: float
  allele2_ru: float

  allele1_reads: int
  allele2_reads: int
  total_reads: int
 
type Genotype = object
  event: Event
  repeat: string
  # total allele sizes in repeat units
  allele1: float
  allele2: float
  # and confidence intervals around the allele size esimates
  GQ: float

type Call = object
  chrom: string
  start: uint
  stop: uint
  # and confidence intervals around position and size
  genotype: Genotype
  quality: float
  # ...

proc tostring*(c: Call): string =
  return &"{c.chrom}\t{c.start}\t{c.stop}\t{c.genotype.allele1:.2f}\t{c.genotype.allele2:.2f}"

# Estimate the size of the smaller allele 
# from reads that span the locus
proc spanning_read_est(reads: seq[Support]): Evidence =
  result.class = "spanning reads"
  result.repeat = reads[0].repeat
  var 
    allele1_bp = NaN
    allele2_bp = NaN
    allele1_ru = NaN
    allele2_ru = NaN
    RepeatCounts: CountTable[uint8]
    Indels: CountTable[int8]

  for read in reads:
    if read.SpanningFragmentLength == 0:
      RepeatCounts.inc(read.SpanningReadRepeatCount)
      Indels.inc(int8(read.SpanningReadCigarInsertionLen) - int8(read.SpanningReadCigarDeletionLen))

  if len(RepeatCounts) >= 2:
    var topRepeatCounts = most_frequent(RepeatCounts, 2)
    allele1_ru = float(topRepeatCounts[0])
    allele2_ru = float(topRepeatCounts[1])
  elif len(RepeatCounts) == 1:
    allele1_ru = float(RepeatCounts.largest[0])
  result.allele1_ru = allele1_ru
  result.allele2_ru = allele2_ru

  if len(Indels) >= 2:
    var topIndels = most_frequent(Indels, 2)
    allele1_bp = float(topIndels[0])
    allele2_bp = float(topIndels[1])
  elif len(Indels) == 1:
    allele1_bp = float(Indels.largest[0])
  result.allele1_bp = allele1_bp
  result.allele2_bp = allele2_bp

# Use a linear model to estimate allele size in bp from sum
# of counts of str repeat units in the anchored reads
# result is in bp insertion from the reference
proc anchored_est(sum_str_counts: uint, depth: float): float =
  #XXX These estimates are from the HTT simulation linear model, need to generalize
  var cofficient = 1.106
  var intercept = 3.348
  var y = log2(float(sum_str_counts)/depth + 1) * cofficient + intercept
  return pow(2,y)

proc genotype*(b:Bounds, tandems: seq[tread], spanners: seq[Support],
              targets: seq[Target], depth: float): Call =
  result.chrom = get_chrom(b.tid, targets)
  result.start = b.left
  result.stop = b.right
  var RUlen = len(b.repeat)

  # Check for spanning reads - indicates a short allele
  if spanners.len == 0:
    result.genotype.allele1 = NaN
  else:
    var small_allele_bp = spanning_read_est(spanners)
    result.genotype.allele1 = small_allele_bp.allele1_bp/float(RUlen)

  # Use anchored reads to estimate long allele
  var sum_str_counts: uint
  for tread in tandems:
    sum_str_counts += tread.repeat_count

  var large_allele_bp = anchored_est(sum_str_counts, depth)
  result.genotype.allele2 = large_allele_bp/float(RUlen)

  # Revise allele esimates using spanning pairs

  # Revise allele esimates using unplaced reads
  # (probably can't be done until all loci are genotyped)

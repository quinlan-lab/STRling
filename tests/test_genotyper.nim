import unittest
include strpkg/genotyper

suite "genotyper suite":

  test "estimate allele size from spanning reads":

    var reads = @[
    Support(SpanningFragmentLength: 0, SpanningReadRepeatCount: 10, SpanningReadCigarInsertionLen: 0, SpanningReadCigarDeletionLen: 0, repeat: "AT", Type: SupportType.SpanningRead),
    Support(SpanningFragmentLength: 0, SpanningReadRepeatCount: 10, SpanningReadCigarInsertionLen: 0, SpanningReadCigarDeletionLen: 0, repeat: "AT", Type: SupportType.SpanningRead),
    Support(SpanningFragmentLength: 0, SpanningReadRepeatCount: 10, SpanningReadCigarInsertionLen: 0, SpanningReadCigarDeletionLen: 0, repeat: "AT", Type: SupportType.SpanningRead),
    Support(SpanningFragmentLength: 0, SpanningReadRepeatCount: 9, SpanningReadCigarInsertionLen: 0, SpanningReadCigarDeletionLen: 2, repeat: "AT", Type: SupportType.SpanningRead),
    ]

    var est = spanning_read_est(reads)
    check est.repeat == "AT"
    check est.allele1_bp == 0.0
    check est.allele2_bp == -2.0
    check est.allele1_ru == 10.0
    check est.allele2_ru == 9.0

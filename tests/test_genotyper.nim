import unittest
include strpkg/genotyper

suite "genotyper suite":

  test "estimate allele size from spanning reads":

    var reads = @[
    Support(SpanningFragmentLength: 0, SpanningReadRepeatCount: 10, SpanningReadCigarInsertionLen: 0, SpanningReadCigarDeletionLen: 0, repeat: "AT"),
    Support(SpanningFragmentLength: 0, SpanningReadRepeatCount: 10, SpanningReadCigarInsertionLen: 0, SpanningReadCigarDeletionLen: 0, repeat: "AT"),
    Support(SpanningFragmentLength: 0, SpanningReadRepeatCount: 10, SpanningReadCigarInsertionLen: 0, SpanningReadCigarDeletionLen: 0, repeat: "AT"),
    Support(SpanningFragmentLength: 0, SpanningReadRepeatCount: 9, SpanningReadCigarInsertionLen: 0, SpanningReadCigarDeletionLen: 2, repeat: "AT"),
    ]

    check spanning_read_est(reads) == Evidence(class: "spanning reads", repeat: "AT", allele1_bp: 0.0, allele2_bp: -2.0, allele1_ru: 10.0, allele2_ru: 9.0)

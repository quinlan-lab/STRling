// Load system configuration and other settings
//load 'pipeline-config.groovy'

// Load Bpipe pipeline stages
load 'pipeline-stages.groovy'

// Define all chromosomes you wish to analyze. They will be prepended by "chr". See https://docs.bpipe.org/Language/Chr/
all_chroms = chr(1..22, 'X','Y')

run {
    str_index +
    "%.${input_type}" * [str_extract] +
    all_chroms * [str_merge_chrom] + str_merge_collect +
    "%.bin" * [str_call_joint] + outliers
}

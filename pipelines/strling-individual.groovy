// Load system configuration and other settings
//load 'pipeline-config.groovy'

// Load Bpipe pipeline stages
load 'pipeline-stages.groovy'

run {
    str_index +
    "%.${input_type}" * [str_extract + str_call_individual]
}

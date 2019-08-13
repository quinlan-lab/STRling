// Bpipe pipeline to simulate paired end reads from a fasta file, with no stutter
// Simulates reads from STR alleles at a provided set of loci

// Set up to run on the redwood cluster

load "pipeline_config.groovy"
load "sim_shared.groovy"

/////////////////////////////
// Run pipeline

run {
        generate_alleles +

        "%.fasta" * [
            generate_reads +
            align_bwa + index_bam + str
        ]
}

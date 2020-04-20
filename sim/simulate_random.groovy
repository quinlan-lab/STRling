// Bpipe pipeline to simulate paired end reads from a fasta file, with no stutter
// Simulates reads from STR alleles at a provided set of loci

load "pipeline_config.groovy"
load "sim_shared.groovy"

mutate_locus = {
    doc """Generate a random heterozygous coding STR loci in the potentially
        pathogenic range."""

    output.dir = "sim"

    produce('*.bed') {
        exec """
            $PYTHON $TOOLS/random_str_alleles.py
                --locus "7 117143769 117143769 CAG"
                --out $dir/
                --num 300
                --min 0
                --max 600
                --fixed 0
                --seed 7
        """
    }
}

// hist based on:
// /uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/data/STR_true_pos/D09-903.cram
@transform("bam")
simulate_str_reads = {
    doc """ """

    output.dir = "sim"

    exec """
        ~/storage/git/STRling/src/strpkg/simulate_reads
            --fasta $REF
            --output $output.prefix
            $input.cram
            $input.bed
    """
}

combine = {

        exec """
            $PYTHON $TOOLS/combine_random_sim_results.py --bed_dir sim --str_dir str --out sim
        """
}

/////////////////////////////
// Run pipeline

run {
        mutate_locus +
        "%.bed" * [
            simulate_str_reads +
            str_extract
        ] + str_merge + 
        "%.bin" * [ str_call ] +
        combine
}

// Bpipe pipeline to simulate paired end reads from a fasta file, with no stutter
// Simulates reads from STR alleles at a provided set of loci

load "pipeline_config.groovy"
load "sim_shared.groovy"

mutate_locus = {
    doc """Generate a random heterozygous coding STR loci in the potentially
        pathogenic range."""

    produce('HTT_rand.bed') {
        exec """
            $PYTHON $TOOLS/random_str_alleles.py
                --locus "4   3076604 3076695 CAG"
                --out $output.bed
                --num 1
                --min -5
                --max 200
                --fixed 0
                --seed 7
        """
    }
}

@transform("bam")
simulate_str_reads = {
    doc """ """

    output.dir = "sim"

    exec """
        /uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/git/str-dev/src/strpkg/simulate_reads
            --fasta $REF
            --output $output.prefix
            /uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/data/STR_true_pos/D09-903.cram
            $input.bed
    """
}

combine = {

        exec """
            $PYTHON $TOOLS/combine_random_sim_results.py --bed_dir sim --str_dir str --out HTT
        """
}

/////////////////////////////
// Run pipeline

run {
        mutate_locus +
        "%.bed" * [
            simulate_str_reads +
            str
        ] + combine
}

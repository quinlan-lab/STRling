// Bpipe pipeline to simulate paired end reads from a fasta file, with no stutter
// Simulates reads from STR alleles at a provided set of loci

// Adjust simulation parameters
PLATFORM="illumina"
total_coverage = 30

def get_fname(path) {
    def x = path.split("/")[-1]
    return(x)
}

/////////////////////////////
// Produce mutated fasta files
generate_alleles = {
    doc "Generate fasta files of STR variants based on the input bed file"

    output.dir = "sim"

    def bedname = get_fname(input.bed)

        produce(bedname.prefix + ".truth.vcf", "*.fasta") {

            exec """
                $PYTHON $TOOLS/generate_str_alleles.py $REF $input.bed --output $dir/ --truth $output.vcf --flank 50000 --id
        """
        }
}

/////////////////////////////
// Generate reads
//generate_reads = {
//    doc "Sample reads from the altered reference sequence segment"
//
//    output.dir = "sim"
//
//    produce( get_fname(input.fasta.prefix) + "_L001_R1.fq", get_fname(input.fasta.prefix) + "_L001_R2.fq") {
//
//        // Set target coverage
//        def coverage = total_coverage
//        def outname = output.prefix[0..-2]
//        exec """
//            $ART/art_illumina -i $input.fasta -p -na
//                -l 150 -ss HS25 -f $coverage
//                -m 350 -s 120
//                -o $outname
//        """
//    }
//}

/////////////////////////////
// Generate reads
generate_reads = {
    doc "Sample reads from the altered reference sequence segment using an empirical insert size distribution"

    output.dir = "sim"

    produce( get_fname(input.fasta.prefix) + "_L001_read1.fq.gz", get_fname(input.fasta.prefix) + "_L001_read2.fq.gz") {

        // Set target coverage
        def coverage = total_coverage
        def outname = output.prefix.prefix[0..-7]
        exec """
            $PYTHON27  ~/storage/git/neat-genreads/genReads.py
                -R 150
                --pe-model $MODELS/fraglen.p
                --gz
                -c $coverage
                -r $input.fasta
                -o $outname

        """
    }
}

/////////////////////////////
// Align reads
@preserve("*.bam")
align_bwa = {
    doc "Concatenate with background reads then align with bwa mem algorithm."

    output.dir = "sim"

    def fname = get_fname(input1)
    def lane = "001"
    def sample = branch.name
    from("fq.gz", "fq.gz") produce(fname.prefix.prefix + ".bam") {
        exec """
            bwa mem -M
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample}\\tSM:${sample}"
            $REF
            $input1.gz
            $input2.gz |
            samtools view -bSuh - | samtools sort -o $output.bam -T $output.bam.prefix
        """, "bwamem"
    }
}

index_bam = {

    output.dir = "sim"

    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

str_extract = {

    output.dir = "str"

    def bamname = get_fname(input.bam)
    def str_ref = get_fname(REF) + ".str"
    produce(bamname.prefix + ".str.bin") {
        exec """
            $STR_NIM extract
                -f $REF
                -g $str_ref
                -v
                $input.bam
                $output.bin
        """
    }
}

str_merge = {

    output.dir = "str"

    from("*.bin") produce("strling-bounds.txt") {
        exec """
            $STR_NIM merge
                -f $REF
                -v
                -o ${output.dir + '/strling'}
                $inputs.bin
        """
    }
}

str_call = {

    output.dir = "str"

    def sample = branch.name.prefix

    from (sample + ".str.bin", sample + ".bam", "strling-bounds.txt") produce(
            sample + "-reads.txt", sample + "-bounds.txt",
            sample + "-spanning.txt", sample + "-unplaced.txt",
            sample + "-genotype.txt") {
    
        def bamname = get_fname(input.bam)
        exec """
            $STR_NIM call
                -v
                -b $input.txt
                -o ${output.dir + '/' + bamname.prefix}
                $input.bam
                $input.bin
        """
    }
}

combine = {

        exec """
            python $TOOLS/combine_random_sim_results.py --bed_dir sim --str_dir str --out HTT
        """
}



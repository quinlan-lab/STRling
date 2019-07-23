// Bpipe pipeline to simulate paired end reads from a fasta file, with no stutter
// Simulates reads from STR alleles at a provided set of loci

// Set up to run on the redwood cluster

INSTALLDIR="~/storage/git/str-dev"
ART="~/tools/art_bin_MountRainier"
REF="/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/ref-data/GATK_Bundle/human_g1k_v37_decoy.fasta"
TOOLS="$INSTALLDIR/sim" //custom R/python scripts
GATK="~/tools/gatk-4.1.2.0/gatk"

//Compiled in regular mode
//STR_NIM="/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/git/str-dev/str"
//Compiled in debug mode
STR_NIM="/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/git/str-dev/src/str"
//Quick share debug mode
//STR_NIM="/tmp/str"

DECOY_REF="$INSTALLDIR/simulate_reads/reference-data/hg19.STRdecoys.sorted.fasta"
EXOME_TARGET="/group/bioi1/harrietd/ref-data/hg19_RefSeq_coding.sorted.bed"

// Adjust simulation parameters
PLATFORM="illumina"
total_coverage = 30

def get_fname(path) {
    def x = path.split("/")[-1]
    return(x)
}

//@preserve("*.bed")
//mutate_locus = {
//    doc """Generate a random heterozygous coding STR loci in the potentially
//        pathogenic range."""
//
//    output.dir = "sim_bed"
//    branch.simID = branch.name
//
//    produce(branch.simID + ".bed") {
//        exec """
//            $INSTALLDIR/simulate_reads/STR_simulation_script.R
//                -L $LOCUS
//                $INSTALLDIR/simulate_reads/reference-data/hg19.simpleRepeat.txt.gz
//                $INSTALLDIR/simulate_reads/reference-data/str-stats
//                -O $output.bed
//                -m 500
//        """
//
//    }
//}

@filter("sorted")
sort_bed = {
    doc "sort bed file"
    output.dir = "sim_bed"
    branch.source_bed = input.bed

    preserve("*.bed") {
        exec """
            bedtools sort -i $input.bed > $output.bed
        """
    }
}

/////////////////////////////
// Produce mutated fasta files
generate_alleles = {
    doc "Generate a VCF of STR mutations and stutter, along with their probabilities"
    output.dir = "vcf_bed"

    def bedname = get_fname(input.bed)

        produce(bedname.prefix + ".truth.vcf", "*.fasta") {

            exec """
                python $TOOLS/generate_str_alleles.py $REF $input.bed --truth $output.vcf --output ${output.prefix.prefix}. --flank 10000
        """
        }
}

/////////////////////////////
// Generate reads
generate_reads = {
    doc "Sample reads from the altered reference sequence segment"
    output.dir = "fastq"

    produce( get_fname(input.fasta.prefix) + "_L001_R1.fq", get_fname(input.fasta.prefix) + "_L001_R2.fq") {

        // Set target coverage
        def coverage = total_coverage
        def outname = output.prefix[0..-2]
        exec """
            $ART/art_illumina -i $input.fasta -p -na
                -l 150 -ss HS25 -f $coverage
                -m 500 -s 50
                -o $outname
        """
    }
}

gzip = {
    output.dir = "fastq"

    produce(input.prefix + ".fastq.gz") {
        exec "cat $inputs.fq | gzip -c > $output.gz"
    }
}


/////////////////////////////
// Align reads
@preserve("*.bam")
align_bwa = {
    doc "Concatenate with background reads then align with bwa mem algorithm."

    def fname = get_fname(input1)
    def lane = "001"
    def sample = branch.name
    from("fastq.gz", "fastq.gz") produce(fname.prefix.prefix + ".bam") {
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
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

str = {

    produce(input.prefix + "-reads.txt", input.prefix + "-bounds.txt") {
        exec """
            $STR_NIM 
                -p 0.7
                -v
                -o $input.prefix
                $input.bam
        """
    }
}

/////////////////////////////
// Run pipeline

run {
        //sort_bed +
        generate_alleles +

        "%.fasta" * [
            generate_reads
        ] +

        "%.fq" * [
            gzip
        ] +

        "%_R*.fastq.gz" * [
            align_bwa + index_bam + str
        ]
}

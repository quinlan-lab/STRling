STRLING='strling'

if(args.any { it.endsWith('.cram') })
    input_type = 'cram'
else
    input_type='bam'

def get_fname(path) {
    def x = path.split("/")[-1]
    return(x)
}
 
str_extract = {
    def sample = branch.name
    def str_ref = get_fname(REF) + ".str"
    produce(sample + ".str.bin") {
        exec """
            $STRLING extract
                -f $REF
                -g $str_ref
                -v
                ${input[input_type]}
                $output.bin
        ""","strling"
    }
}

str_merge = {
    from("*.bin") produce("strling-bounds.txt") {
        exec """
            $STRLING merge
                -f $REF
                -v
                $inputs.bin
        ""","strling"
    }
}

str_call_individual = {
    def sample = branch.name.prefix
    from (sample + ".str.bin", sample + "." + input_type) produce(
            sample + "-bounds.txt", sample + "-unplaced.txt",
            sample + "-genotype.txt") {
        exec """
            $STRLING call
                -f $REF
                -v
                -o $sample
                ${input[input_type]}
                $input.bin
        ""","strling"
    }
}

str_call_joint = {
    def sample = branch.name.prefix
    from (sample + ".str.bin", sample + "." + input_type, "strling-bounds.txt") produce(
            sample + "-bounds.txt", sample + "-unplaced.txt",
            sample + "-genotype.txt") {
        exec """
            $STRLING call
                -f $REF
                -v
                -b $input.txt
                -o $sample
                ${input[input_type]}
                $input.bin
        ""","strling"
    }
}


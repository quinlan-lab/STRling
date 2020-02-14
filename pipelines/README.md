# Bpipe pipelines

Place the strling executable in your PATH, or set `-p STRLING=/path/to/strling`

## Individual calling
`bg-bpipe run -p REF=/path/to/ref-genome.fasta strling-individual.groovy *.cram`

## Joint calling
`bg-bpipe run -p REF=/path/to/ref-genome.fasta strling-joint.groovy *.cram`

Note that nim does not play nicely with bpipe's ctl-C handing, so have to use the background bpipe command `bg-bpipe` instead of the regular `bpipe`.

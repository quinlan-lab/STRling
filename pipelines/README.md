# Bpipe pipelines

Note that nim does not play nicely with bpipe's ctl-C handing, so have to use the background bpipe command `bg-bpipe` instead of the regular `bpipe`.

Place the strling executable in your PATH, or set `-p STRLING=/path/to/strling`

`bg-bpipe run strling-joint.groovy -p REF=/path/to/ref-genome.fasta pipeline.groovy *.cram`

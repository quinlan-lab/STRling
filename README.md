[![Build Status](https://travis-ci.org/quinlan-lab/STRling.svg?branch=master)](https://travis-ci.org/quinlan-lab/STRling)

## Install from source

Install nim:  
`curl https://nim-lang.org/choosenim/init.sh -sSf > init.sh && sh init.sh`

Install STRling:  
```
git clone <URL>
cd STRling
nimble install
```

Compile options for development:  

Compile in fast mode (danger) with read names reported:  
`nim c -d:danger -d:release src/strling.nim`

## Run

#### extract informative pairs to a binary format
```
name=hg002
strling extract -f $reference_fasta /path/to/$name.cram $name.bin
```

#### call strs on the extract binary data

```
mkdir -p str-results/
strling call --output-prefix str-results/$name -f $reference_fasta /path/to/$name.cram $name.bin
```

## Outputs

The main output file is `{$prefix}-genotype.txt`. It reports all STR expansion loci that pass thresholds as well as any provided as input. The columns are:
- chrom: chromosome/contig name
- left: predicted left boundry of STR locus
- right: predicted right boundry of STR locus
- repeatunit: predicted STR repeat unit
- allele1\_est: estimated size of the shorter allele supported by spanning reads (if any)
- allele2\_est: estimated size of the larger allele supported by anchored reads
- total\_reads: number of reads supporting an expansion at this locus
- spanning\_reads: number of reads that span the locus
- spanning\_pairs: number of read pairs that span the locus
- left\_clips: number of soft-clipped reads supporting the left side of the locus position
- right\_clips: number of soft-clipped reads supporting the right side of the locus position
- unplaced\_pairs: number of unplaced STR reads assigned to this locus (will only be >0 for a uniquely expanded repeat unit)
- depth: local median depth around the locus
- sum\_str\_counts: the sum of STR repeat units in all reads supporting an expansion

Some additional outputs are provided with detailed supporting evidence used to make the genotype calls:
- Putative str bounds: `{$prefix}-bounds.txt`
- All str-like reads: `{$prefix}-reads.txt`
- Spanning reads and spanning pairs:`{$prefix}-spanning.txt`
- Counts of str-like reads that are unplaced (could not be assigned to a locus): '{$prefix}-unplaced.txt`


## Run tests
`nimble tests`

If you get the error:
`could not load: libhts.so`

try:  
`export LD_LIBRARY_PATH=./htslib/`



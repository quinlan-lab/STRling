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
`nim c -d:danger -d:qname src/strling.nim`

## Run

#### extract informative pairs to a binary format
```
name=hg002
strling extract -v -f $reference_fasta /path/to/$name.cram $name.bin
```

#### call strs on the extract binary data

```
mkdir -p str-results/
strling call --output-prefix str-results/$name -f $reference_fasta /path/to/$name.cram $name.bin
```


## Run tests
`nimble tests`

If you get the error:
`could not load: libhts.so`

try:  
`export LD_LIBRARY_PATH=./htslib/`


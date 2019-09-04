[![Build Status](https://travis-ci.com/hdashnow/str-dev.svg?branch=master)](https://travis-ci.com/hdashnow/str-dev)

## Install from source

Install nim:  
`curl https://nim-lang.org/choosenim/init.sh -sSf > init.sh && sh init.sh`

Install str:  
```
git clone <URL>
cd str-dev
nimble install
```

Compile options for development:  

Compile in fast mode (danger) with read names reported:  
`nim c -d:danger -d:qname src/str.nim`

## Run

#### extract informative pairs to a binary format
```
name=hg002
str extract -v -f $reference_fasta /path/to/$name.cram $name.bin
```

#### call strs on the extract binary data

```
mkdir -p str-results/
str call --output-prefix str-results/$name -f $reference_fasta /path/to/$name.cram $name.bin
```


## Run tests
`nimble tests`

If you get the error:
`could not load: libhts.so`

try:  
`export LD_LIBRARY_PATH=./htslib/`


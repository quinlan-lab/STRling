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
`./str`

## Run tests
`nimble tests`

If you get the error:
`could not load: libhts.so`

try:  
`export LD_LIBRARY_PATH=./htslib/`


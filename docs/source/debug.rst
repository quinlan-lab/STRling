Debug
=====

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Common issues
-------------

- If you get the error:
`could not load: libhts.so`

Try:
`export LD_LIBRARY_PATH=./htslib/`

- Empty -bounds.txt/-genotype.txt:

Check the mapping quality (MAPQ) of your reads.
STRling requires a default minimum MAPQ of 40 for non-repetitive reads.
If the MAPQs are low, and you'd still like to consider these reads,
you can change the MAPQ setting using: 
`-q, --min-mapq=MIN_MAPQ`  
Ensure to do this at all workflow stages where it is an option.


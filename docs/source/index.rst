.. STRling documentation master file, created by
   sphinx-quickstart on Sun Dec 22 16:03:09 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to STRling's documentation!
===================================

STRling code: https://github.com/quinlan-lab/STRling

STRling (pronounced like “sterling”) is a method to detect large STR expansions from short-read sequencing data. It is capable of detecting novel STR expansions, that is expansions where there is no STR in the reference genome at that position (or a different repeat unit from what is in the reference). It can also detect STR expansions that are annotated in the reference genome. STRling uses kmer counting to recover mis-mapped STR reads. It then uses soft-clipped reads to precisely discover the position of the STR expansion in the reference genome.

Dashnow, H., Pedersen, B.S., Hiatt, L. et al. STRling: a k-mer counting approach that detects short tandem repeat expansions at known and novel loci. Genome Biol 23, 257 (2022). https://doi.org/10.1186/s13059-022-02826-4

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   run
   outputs
   workflows
   test
   debug
   contribute
   cite


Quick start
===========


Install
-------

We recommending downloading the static binary.

Download the `strling` binary from the latest release `here <https://github.com/quinlan-lab/STRling/releases/latest>`_.

Make it executable:
`chmod +x strling`

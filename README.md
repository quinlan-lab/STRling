![STRling logo](docs/strling-logo-webres.png)

[![CI](https://github.com/quinlan-lab/STRling/actions/workflows/ci.yml/badge.svg)](https://github.com/quinlan-lab/STRling/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/strling/badge/?version=latest)](https://strling.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

__[STRling Documentation](https://strling.readthedocs.io/en/latest/)__

STRling (pronounced like “sterling”) is a method to detect large STR expansions from short-read sequencing data. It is capable of detecting novel STR expansions, that is expansions where there is no STR in the reference genome at that position (or a different repeat unit from what is in the reference). It can also detect STR expansions that are annotated in the reference genome. STRling uses kmer counting to recover mis-mapped STR reads. It then uses soft-clipped reads to precisely discover the position of the STR expansion in the reference genome.

## Install and Run STRling

Please see the [STRling Documentation](https://strling.readthedocs.io/en/latest/) for installation and running instructions.

## Citation

For more details able the algorithm check out our paper.

If using STRling, please cite: 

Dashnow, H., Pedersen, B.S., Hiatt, L. et al. STRling: a k-mer counting approach that detects short tandem repeat expansions at known and novel loci. Genome Biol 23, 257 (2022). <https://doi.org/10.1186/s13059-022-02826-4>

[![Build Status](https://travis-ci.org/quinlan-lab/STRling.svg?branch=master)](https://travis-ci.org/quinlan-lab/STRling)
[![Documentation Status](https://readthedocs.org/projects/strling/badge/?version=latest)](https://strling.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

__[STRling Documentation](https://strling.readthedocs.io/en/latest/)__

STRling is still in development. Please report bugs via GitHub issues.

STRling (pronounced like “sterling”) is a method to detect large STR expansions from short-read sequencing data. It is capable of detecting novel STR expansions, that is expansions where there is no STR in the reference genome at that position (or a different repeat unit from what is in the reference). It can also detect STR expansions that are annotated in the reference genome. STRling uses kmer counting to recover mis-mapped STR reads. It then uses soft-clipped reads to precisely discover the position of the STR expansion in the reference genome.

## Install and Run STRling

Download the latest release of `strling` from [here](https://github.com/quinlan-lab/STRling/releases/latest).

Please see the [STRling Documentation](https://strling.readthedocs.io/en/latest/) for running instructions.

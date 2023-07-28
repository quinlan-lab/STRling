Filter
=====

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Filtering output
----------------

In our experience, homopolymers, and some repetitive regions are enriched for false positives. You may find the following filters useful:

Exclude homopolymers
Require the repeatunit column be greater than 1 character.
STRs.tsv files: `awk '(length($6) > 1)'`
*-genotype.txt files: `awk '(length($4 > 1)'`

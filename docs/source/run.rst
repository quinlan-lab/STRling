Run
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Creating a genome STR index (optional)
--------------------------------------

Creates a bed file of large STR regions in the reference genome. This step is peformed automatically as part of `strling extract`. However, when running multiple samples, it is more efficient to do it once, then pass the file to `strling extract` using the `-g` option.

.. code-block:: bash

    strling index $reference_fasta

Single sample
-------------

Extract informative pairs to a binary format.

.. code-block:: bash

  strling extract -f $reference_fasta /path/to/$name.cram $name.bin

Call strs on the extract binary data.

.. code-block:: bash

  mkdir -p str-results/
  strling call --output-prefix str-results/$name -f $reference_fasta /path/to/$name.cram $name.bin


Joint calling
-------------

Extract informative read pairs to a binary format for a single sample (same as above, you can use the same bin files).

.. code-block:: bash

  strling extract -f $reference_fasta /path/to/$sample2.cram $sample1.bin
  strling extract -f $reference_fasta /path/to/$sample2.cram $sample2.bin

Joint call str loci across all samples. Requires minimum read evidence from at least one sample.

.. code-block:: bash

  mkdir -p str-results/
  strling merge --output-prefix str-results/joint -f $reference_fasta $sample1.bin $sample2.bin

Call genotypes/estimate allele sizes for all loci in each sample.

.. code-block:: bash

  strling call --output-prefix str-results/$sample1 -b str-results/joint-bounds.txt -f $reference_fasta /path/to/$sample1.cram $sample1.bin
  strling call --output-prefix str-results/$sample2 -b str-results/joint-bounds.txt -f $reference_fasta /path/to/$sample2.cram $sample2.bin

Find outliers: loci that are expanded in one individual relative to other individuals in the joint called cohort.

.. code-block:: bash

  python scripts/strling-outliers.py --genotypes *-genotype.txt --unplaced *-unplaced.txt

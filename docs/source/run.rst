Run
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Creating a genome STR index (optional)
--------------------------------------

Creates a bed file of large STR regions in the reference genome. This step is performed automatically as part of `strling extract`. However, when running multiple samples, it is more efficient to do it once, then pass the file to `strling extract` using the `-g` option.

.. code-block:: bash

    strling index $reference_fasta

Single sample
-------------

**Extract informative pairs to a binary format**

.. code-block:: bash

  strling extract -f $reference_fasta /path/to/$sample.cram $sample.bin

Output file: $sample.bin - a binary file describing STR-containing reads

**Call strs on the extract binary data**

.. code-block:: bash

  mkdir -p str-results/
  strling call --output-prefix str-results/$sample -f $reference_fasta /path/to/$sample.cram $sample.bin

Output files:
$sample-bounds - STR loci interrogated in that sample
$sample-genotype.txt - Locus size estimates and other per-locus information (see :ref:`outputs`)
$sample-unplaced.txt - Counts of unplaced STR reads that could not be assigned to a specific locus

Joint calling
-------------

**Extract informative read pairs to a binary format for each sample**

This step is the same as above, you can use the same bin files.

.. code-block:: bash

  strling extract -f $reference_fasta /path/to/$sample2.cram $sample1.bin
  strling extract -f $reference_fasta /path/to/$sample2.cram $sample2.bin

Output file(s): $sample1.bin, $sample2.bin, ... - binary files describing STR-containing reads

**Joint call str loci across all samples**

Requires minimum read evidence from at least one sample.

.. code-block:: bash

  mkdir -p str-results/
  strling merge --output-prefix str-results/joint -f $reference_fasta $sample1.bin $sample2.bin

Output file: joint-bounds.txt - positions of STR loci found by combining across all individuals, used for the call stage when joint calling

Peak memory usage for joint calling increases linearly by ~63 MB per 30-40X human WGS. This can be prohibitively high for large cohorts. To reduce memory requirements, analysis can be parallelized by chromosome using `--chromosome` then subsequently merged. See the workflows for examples. This reduces memory usage to ~5 MB/sample.

**Call genotypes/estimate allele sizes for all loci in each sample**

.. code-block:: bash

  strling call --output-prefix str-results/$sample1 -b str-results/joint-bounds.txt -f $reference_fasta /path/to/$sample1.cram $sample1.bin
  strling call --output-prefix str-results/$sample2 -b str-results/joint-bounds.txt -f $reference_fasta /path/to/$sample2.cram $sample2.bin

Output files as above.

**Find outliers**

Finds loci that are expanded in one individual relative to other individuals in the joint called cohort.

.. code-block:: bash

  strling-outliers.py --genotypes *-genotype.txt --unplaced *-unplaced.txt

Output files:
STRs.tsv - a single file with all loci in all samples and their outlier p-values (see :ref:`outputs`)
$sample.STRs.tsv - the same data, filtered to a single individual

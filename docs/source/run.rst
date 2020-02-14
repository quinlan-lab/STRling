Run
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Single sample
-------------

Extract informative pairs to a binary format

.. code-block:: bash

  strling extract -f $reference_fasta /path/to/$name.cram $name.bin

Call strs on the extract binary data

.. code-block:: bash

  mkdir -p str-results/
  strling call --output-prefix str-results/$name -f $reference_fasta /path/to/$name.cram $name.bin


Joint calling
-------------

Extract informative pairs to a binary format for a single sample (same as above, you can use the same bin files)

.. code-block:: bash

  strling extract -f $reference_fasta /path/to/$sample2.cram $sample1.bin
  strling extract -f $reference_fasta /path/to/$sample2.cram $sample2.bin

Joint call str loci across all samples

.. code-block:: bash

  mkdir -p str-results/
  strling merge --output-prefix str-results/joint -f $reference_fasta $sample1.bin $sample2.bin

Call estimate allele sizes for each sample

.. code-block:: bash

  strling call --output-prefix str-results/$sample1 -b str-results/joint-bounds.txt -f $reference_fasta /path/to/$sample1.cram $sample1.bin
  strling call --output-prefix str-results/$sample2 -b str-results/joint-bounds.txt -f $reference_fasta /path/to/$sample2.cram $sample2.bin

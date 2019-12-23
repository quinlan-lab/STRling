Run tests
=========

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Run unit tests
--------------

From the `STRling directory` run:

`nimble test`

End-to-end tests
----------------

Download data for end to end tests:

.. code-block:: bash

  mkdir -p run_tests
	cd run_tests
	wget -O data.zip https://ndownloader.figshare.com/articles/11367851?private_link=207b2a78ddcbdb28781b
	unzip data.zip
	rm data.zip
  # Download reference
	wget -O chr4.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr4.fa.gz
	gzip -d chr4.fa.gz
  # End-to-end test single
	mkdir -p single
	../strling extract -f chr4.fa -g chr4.fa.str -v CANVAS.chr4-39348425_AAGGG_0_400.bam single/CANVAS.chr4-39348425_AAGGG_0_400.str.bin
	../strling call -v -o single/CANVAS.chr4-39348425_AAGGG_0_400 CANVAS.chr4-39348425_AAGGG_0_400.bam single/CANVAS.chr4-39348425_AAGGG_0_400.str.bin
  # End-to-end test joint
	mkdir -p joint
	for i in *.bam; do ../strling extract -f chr4.fa -g chr4.fa.str -v $i joint/$i.str.bin; done
	../strling merge -f chr4.fa -v -o joint/strling joint/*.bin
	for i in *.bam; do ../strling call -v -b joint/strling-bounds.txt -o joint/$i $i joint/$i.str.bin; done

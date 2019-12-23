Outputs
=======

.. toctree::
   :maxdepth: 2
   :caption: Contents:


The main output file is `prefix-genotype.txt`. It reports all STR expansion loci that pass thresholds as well as any provided as input. The columns are:

* chrom: chromosome/contig name
* left: predicted left boundry of STR locus
* right: predicted right boundry of STR locus
* repeatunit: predicted STR repeat unit
* allele1\_est: estimated size of the shorter allele in repeat units relative to the reference, from spanning reads (if any). "na" indicats no reads support an allele shorter than the read length, so both may be large.
* allele2\_est: estimated size of the larger allele in repeat units relative to the reference, from anchored reads
* total\_reads: number of reads supporting an expansion at this locus
* spanning\_reads: number of reads that span the locus
* spanning\_pairs: number of read pairs that span the locus
* left\_clips: number of soft-clipped reads supporting the left side of the locus position
* right\_clips: number of soft-clipped reads supporting the right side of the locus position
* unplaced\_pairs: number of unplaced STR reads assigned to this locus (will only be >0 for a uniquely expanded repeat unit)
* depth: local median depth around the locus
* sum\_str\_counts: the sum of STR repeat units in all reads supporting an expansion

Some additional outputs are provided with detailed supporting evidence used to make the genotype calls:

* Putative str bounds: `prefix-bounds.txt`
* Counts of str-like reads that are unplaced (could not be assigned to a locus): `prefix-unplaced.txt`

Only output when compiled with `-d:debug`:

* All str-like reads: `prefix-reads.txt`
* Spanning reads and spanning pairs:`prefix-spanning.txt`

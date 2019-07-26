#!/usr/bin/env python
"""Output:
A single vcf file containing the true genotype that is being simulated.
A bed file corresponding to the region around the STR.
"""

import sys
from argparse import (ArgumentParser, FileType)
from collections import OrderedDict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import pysam 
import random
import pandas as pd
import itertools
__author__ = 'Harriet Dashnow'
__credits__ = ['Harriet Dashnow']
__license__ = 'MIT'
__version__ = '0.1.0'
__email__ = 'h.dashnow@gmail.com'

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = ArgumentParser(description=('Produce fasta files and VCFs'
    ' frequencies for a given set of STR loci and alleles. Also provides a bed file for each'
    ' locus defining a region around that locus.'))
    parser.add_argument(
        'ref', type=str,
        help='Fasta reference')
    parser.add_argument(
        'bed', type=str,
        help='bed file containing genomic locations of STRs and their repeat units. Genomic locations should be relative to the fasta reference. format: chr start stop name, where name is in the format repeatunit_genotype, e.g. CAG_-2/1')
    parser.add_argument(
        '--output', type=str, required=False, default='',
        help='Base name for output files, including vcfs and bed files.')
    parser.add_argument(
        '--id', action="store_true",
        help='Prefix individual fasta and bed output files with a numerical id.')
    parser.add_argument(
        '--truth', type=str, required=False, default='truth.vcf',
        help='File name for output vcf of true genotypes for all loci. (default: %(default)s)')
    parser.add_argument(
        '--flank', type=int, default=10000,
        help='Number of flanking bases to include in the output on either side of the STR. (default: %(default)s)')
    parser.add_argument(
        '--target', type=str,
        help='bed file containing genomic locations of the region to the simulated. Warning: variants outside these regions will be excluded.')
    parser.add_argument(
        '--seed', required=False,
        help='Random seed (can be any hashable input).')
    return parser.parse_args()


def circular_permuted(x):
    """Generates all possible circular permutations of the input.

    Args:
        x (str or any iterable?)

    Returns:
        list: All circular permutations of x
    """
    return([x[i:] + x[:i] for i in range(len(x))])

def self_and_rev_complement(in_dna):
    """Returns the input DNA sequence and its reverse complement

    Args:
        in_dna (str): valid DNA bases e.g. ACTGN

    Returns:
        list: [in_dna, reverse_complement]
    """
    all_possible = [in_dna]
    # Get reverse complement
    dna = Seq(in_dna, generic_dna)
    rev_complement = str(dna.reverse_complement())
    all_possible.append(rev_complement)
    return(all_possible)

def normalise_str(in_dna):
    """Find all possible eqivalent Short Tandem Repeat (STR) DNA sequences.
    And return the first alphabetically.
    For example, TA = AT. But would return AT.

    Args:
        in_dna (str): STR repeat unit consisting of valid DNA bases e.g. ACTGN

    Returns:
        str: normalised STR sequence
    """
    all_possible = []
    # Circularly permute original sequence and reverse complement
    for seq in self_and_rev_complement(in_dna):
        for permuted_seq in circular_permuted(seq): # Switch to faster permutation (6)
            all_possible.append(permuted_seq)

    # Sort and take the first
    all_possible.sort()
    return(all_possible[0])

def is_dna(a):
    """Return True if input contains only DNA characters, otherwise return False

    Args:
        a (str): string to test

    Returns:
        bool: True if all characters are DNA, otherwise False
    """
    if len(a) == 0:
        return(False)
    dna_chars = 'atcgnATCGN'
    return all(i in dna_chars for i in a)

def parse_bed(bedfilename, position_base = 0, bed_dict = OrderedDict(), pad_left = 1):
    """Parse regions from bed file. Ignore lines starting with #.

    Args:
        bedfilename (str): bed format text file in in
        format: chr start stop name [...] (additional columns ignored)
        where name is in the format repeatunit_genotype, e.g.
        CAG_-2/1 (CAG repeat unit, 2 CAG deleation/1 CAG insertion)
        AT_0/3 (AT repeat unit, same as reference/3 AT insertion)
        bed_dict (dict): existing genomic regions to include in the output.
            format: bed_dict[unique_id] = {'chr':str, 'start':int, 'stop':int, 'name':None}
        position_base (int): 0 or 1. The starting position the genomic regions
            are measured relative to in the input file. i.e. the numbering of the
            first base in the reference genome. If 1, all postions will
            be converted to base-0. Assumed to be 0 by default.
        pad_left (int): Minus pad_left to the start of each locus. This provides
            this many bases of padding to the start(left) of each locus so that
            all repeat units can be deleated without leaving a blank alt allele.
            Setting this to 0 could result in blank alt alleles, and therefore
            generate an invalid VCF. Default = 1.

    Returns:
        dict: bed_dict[unique_id] = {'chr':str, 'start':int, 'stop':int,
                                    'name':str or None, 'repeatunit':str or None,
                                    'deltas': [int, int] or None}
            Genomic regions, in base-0 (i.e. bed format).
    """
    with open(bedfilename) as bedfile:
        if position_base == 0:
            base_shift = 0
        else:
            base_shift = -1
        for bedfile_line in bedfile:
            if bedfile_line.startswith('#') or bedfile_line == '\n':
                continue
            split_line = bedfile_line.split()
            ref_chr = split_line[0]
            ref_start = int(split_line[1])
            ref_stop = int(split_line[2])
            # If repeat start and end are the wrong way around, swap them
            if ref_stop < ref_start:
                ref_start = int(split_line[2])
                ref_stop = int(split_line[1])
                sys.stderr.write('Warning, bed start position greater than end position for line:')
                sys.stderr.write(bedfile_line)
            # Change to base-0 if needed. End remains unchanged as per bed format.
            # Apply padding bases.
            ref_start = ref_start + base_shift - pad_left
            if len(split_line) > 3:
                name = split_line[3]
                split_name = name.split('_')
                # Parse out STR repeat unit and target genotype if present (required?)
                if len(split_name) >= 2:
                    repeatunit = split_name[0]
                    if is_dna(repeatunit):
                        repeatunit = repeatunit.upper()
                    genotype = split_name[1].split('/')
                    if len(genotype) == 2:
                        deltas = [int(x) for x in genotype]
            else:
                name = None
                repeatunit = None
                deltas = None
            unique_id = '{0}-{1}-{2}-{3}'.format(ref_chr, ref_start, ref_stop, name)
            if unique_id not in bed_dict:
                bed_dict[unique_id] = {'chr':ref_chr, 'start':ref_start,
                                        'stop':ref_stop, 'name':name,
                                        'repeatunit':repeatunit,
                                        'deltas': deltas}
    return bed_dict

def mutate_str(ref_sequence, repeatunit, delta, random=False):
    """Mutate a DNA sequence containing a microsatellite by inserting or
    deleating repeat units of that microsatellite.

    Args:
        ref_sequence (str): The reference DNA sequence to be mutated.
        repeatunit (str): DNA repeat unit to be inserted/deleted from the
            ref_sequence.
        delta (int): Number of repeat units to add/remove. Positive if creating
            an insertion, negative if a deletion. 0 will return the input
            sequence.
        random (bool): Insert/delete repeat units from a random position.
        If False, inserts/deletes in the left-most position. Note, GATK seems
        to silently fail on vcfs with indels to the right of an imperfect repeat.

    Returns:
        str: The mutated DNA sequence.
    """
    # Check ref_sequence and repeatunit are DNA
    if not is_dna(ref_sequence):
        raise ValueError("{0} is not a valid DNA sequence".format(ref_sequence))
    if not is_dna(repeatunit):
        raise ValueError("{0} is not a valid DNA sequence".format(repeatunit))
    repeatunitlen = len(repeatunit)
    if delta == 0:
        return(ref_sequence)
    if delta > 0: # Insertion
        i_max = len(ref_sequence)-repeatunitlen
        if random:
            # select a position at random, sampling without replacement
            i_range = random.sample(range(i_max), i_max)
        else:
            i_range = range(i_max)
        for i in i_range:
            # Check that there is a repeat unit to the right of this position
            bases_to_right = ref_sequence[i:i+repeatunitlen]
            if normalise_str(bases_to_right) == normalise_str(repeatunit):
                # Insert repeat units at that position
                # (by replicating the reference version of the repeat unit)
                new_sequence = ref_sequence[:i] + bases_to_right * delta + ref_sequence[i:]
                return(new_sequence)
        # Check the last few bases for the repeat unit
        last_bases = ref_sequence[-repeatunitlen:]
        if normalise_str(last_bases) == normalise_str(repeatunit):
            # Insert repeat units at end of sequence
            new_sequence = ref_sequence + last_bases * delta
            return(new_sequence)
        #XXX TODO If no solution was found, insert to the left, or random?
        raise ValueError("The repeat unit {0} was not found in {1}".format(repeatunit, ref_sequence))
    if delta < 0: # Deletion
        deletion_size = repeatunitlen * -delta
        if deletion_size > len(ref_sequence):
            raise ValueError("Deletion of {0} {1} repeat units is larger than the input sequence {2}".format(-delta, repeatunit, ref_sequence))
        i_max = len(ref_sequence) - deletion_size + 1
        if random:
            i_range = random.sample(range(i_max), i_max)
        else:
            i_range = range(i_max)
        for i in i_range:
            bases_to_right = ref_sequence[i:i+deletion_size]
            # Check if the first few bases contain the repeat unit (or a transposition of it)
            for j in range(len(bases_to_right)):
                ref_seg = bases_to_right[j:j+repeatunitlen]
                # Find the version of the repeat unit present in the ref
                if normalise_str(repeatunit) == normalise_str(ref_seg):
                    ref_unit = ref_seg
                    # Check if there are enough repeat units to be deleted
                    if bases_to_right == ref_unit * -delta:
                        # Generate mutated sequence
                        new_sequence = ref_sequence[:i] + ref_sequence[i+deletion_size:]
                        return(new_sequence)
                else:
                    break
        raise ValueError("There were not {0} copies of {1} repeat unit available to be deleted in {2}.".format(-delta, repeatunit, ref_sequence))

def right_trim_indel(ref, alt):
    """Generate the shortest representation of an indel my removing the rightmost bases.
    XXX Currently only working with a single alt, should be generalised for multiple.

    Args:
        ref (str): The version of the sequence in the reference genome.
        alt (str): Alternative version of the sequence caused by an insertion or
        deletion.

    Returns:
        (ref_normalised str, alt_normalised str)
    """
    if min(len(ref), len(alt)) <= 1:
        return ref, alt
    if ref == alt:
        raise ValueError("Ref and alt are the same. ref: {0} alt: {1}".format(ref, alt))
    for i in range(1, min(len(ref), len(alt))):
        if ref[-i] != alt[-i]: # Check if the rightmost bases are different
            if i == 1: # The last bases are identical, so can't be trimmed
                return ref, alt
            else:
                return ref[:-i+1], alt[:-i+1]
    return ref[:-i], alt[:-i]

def left_trim_indel(ref, alt, pos):
    """Generate the shortest representation of an indel my removing the leftmost bases.
    XXX Currently only working with a single alt, should be generalised for multiple.

    Args:
        ref (str): The version of the sequence in the reference genome.
        alt (str): Alternative version of the sequence caused by an insertion or
        deletion.
        pos (int): Start position of the variant

    Returns:
        (ref_normalised str, alt_normalised str, pos_normalised int)
    """
    if min(len(ref), len(alt)) <= 1:
        return ref, alt, pos
    if ref == alt:
        raise ValueError("Ref and alt are the same. ref: {0} alt: {1}".format(ref, alt))
    for i in range(0, min(len(ref), len(alt))):
        new_pos = pos + i
        if ref[i] != alt[i]: # Check if the leftmost bases are different
            if i == 0: # The last bases are identical, so can't be trimmed
                return ref, alt, new_pos
            else:
                return ref[i:], alt[i:], new_pos
    return ref[i:], alt[i:], new_pos

def replace_variant(ref, variant, start, stop=None):
    """Take a string, ref. Insert a string variant that replaces the bases
    in ref from start to stop, inclusive.
    start and stop are 0-based Pythonic coordinates
    if stop is None, the variant will simply be inserted before the start base
    """
    if stop == None:
        stop = start
    assert stop >= start
    assert start > 0 and stop > 0
    assert start <= len(ref)
    assert stop <= len(ref)
    return ref[:start] + variant + ref[stop:]

def get_vcf_writer(vcf_outfile, samples=['SAMPLE'], source='generate_str_alleles.py', ref=''):
    """Write vcf header to file given. Return writer object for that file.

    Args:
        vcf_outfile (str): filename to write the vcf header to
        samples (list of str): names of samples
        source (str): program used to generate vcf (this script by default)
        ref (str): reference genome

    Return:
        file writer object
    """

    vcf_out = open(vcf_outfile, 'w')

    vcf_out.write('##fileformat=VCFv4.1\n')
    vcf_out.write('##source={0}\n'.format(source))
    vcf_out.write('##reference={0}\n'.format(ref))
    vcf_out.write('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat Unit">\n')
    vcf_out.write('##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference Length of Repeat">\n')
    vcf_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    samples_header = '\t'.join(samples)
    vcf_out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}\n'.format(samples_header))

    return(vcf_out)

def get_alt_genotype(ref, alt1, alt2=None):
    """For a given reference allele, and two alternates, return the appropriate
    ALT and GT (genotype) fields for to be written to the VCF

    Args:
        ref (str): reference allele
        alt1 (str): first alternate allele
        alt2 (str): second alternate allele

    Return:
        tuple (vcf_alt, vcf_gt)
        vcf_alt (str):
        vcf_gt (str):
    """
    if alt2 == None:
        alt2 = alt1
    if alt1.upper() == ref.upper() and alt2.upper() == ref.upper():
        vcf_alt = '.'
        vcf_gt = '0/0'
    elif alt1.upper() == alt2.upper():
        vcf_alt = alt1
        vcf_gt = '1/1'
    elif alt1.upper() == ref.upper():
        vcf_alt = alt2
        vcf_gt = '0/1'
    elif alt2.upper() == ref.upper():
        vcf_alt = alt1
        vcf_gt = '0/1'
    else:
        vcf_alt = ','.join([alt1, alt2])
        vcf_gt = '1/2'
    return(vcf_alt, vcf_gt)

def main():
    # Parse command line arguments
    args = parse_args()
    outfile_base = args.output

    vcf_truth = get_vcf_writer(args.truth)

    if args.seed: #XXX required?
        random.seed(args.seed)

    # Add this many bases of padding to the start (left) of each locus so that
    # all repeat units can be deleated without leaving a blank alt allele.
    pad_left = 1

    if args.target:
        target_dict = parse_bed(args.bed, pad_left=pad_left)

    # Parse STR regions that need to be simulated
    bed_dict = parse_bed(args.bed, pad_left=pad_left)
    # get corresponding bit of the fasta file
    fastafile = pysam.Fastafile(args.ref)

    vcf_probs_dict = {}
    for (region, i) in zip(bed_dict, range(len(bed_dict))): 
        chrom = bed_dict[region]['chr']
        start = bed_dict[region]['start'] # These positions are in base-0
        stop = bed_dict[region]['stop']
        name = bed_dict[region]['name']
        repeatunit = bed_dict[region]['repeatunit']
        deltas = bed_dict[region]['deltas']
        # note fetch step requires zero-based indexing c.f. 1-based indexing in vcf and bed dict
        ref_sequence = fastafile.fetch(chrom, start, stop).upper()
        # Repeat length in number of repeat units (esimtated from input, not actually counted)
        ref_RL = (len(ref_sequence) - pad_left) / len(repeatunit)

        bedout_line = '{0}\t{1}\t{2}\t{3}\n'.format(chrom, start - args.flank, stop + args.flank, name)

        try:
            allele1 = mutate_str(ref_sequence, repeatunit, delta = deltas[0])
            allele2 = mutate_str(ref_sequence, repeatunit, delta = deltas[1])
        except ValueError as e:
            sys.stderr.write(region + '\n')
            sys.stderr.write(str(e) + '\n')
            continue

        # Get fasta with flanking sequence
        this_ref = fastafile.fetch(chrom, start - args.flank, stop + args.flank).upper()
        rel_start = args.flank
        rel_stop = rel_start + stop - start
        # Insert new alleles into the reference fasta with flanks
        allele1_fasta = replace_variant(this_ref, allele1, rel_start, rel_stop)
        allele2_fasta = replace_variant(this_ref, allele2, rel_start, rel_stop)

        if args.id:
            out_prefix = '{}{}.{}-{}_{}_{}_{}'.format(args.output, i, chrom, start, repeatunit, deltas[0], deltas[1])
        else:
            out_prefix = '{}{}-{}_{}_{}_{}'.format(args.output, chrom, start, repeatunit, deltas[0], deltas[1])

        # Write new alleles to fasta
        sequences = [SeqRecord(Seq(allele1_fasta), id='{}_{}'.format(repeatunit, deltas[0]), description=''),
                    SeqRecord(Seq(allele2_fasta), id='{}_{}'.format(repeatunit, deltas[1]), description='')]
        with open(out_prefix + '.fasta', "w") as outfasta:
            SeqIO.write(sequences, outfasta, "fasta")

        # Write bed file
        with open(out_prefix + '.bed', "w") as bed_out:
            bed_out.write(bedout_line)

        # Write the true alleles
        vcf_start = start + 1 # convert to base 1 for vcf file
        vcf_alt, vcf_gt = get_alt_genotype(ref_sequence, allele1, allele2)
        vcf_id = '.'
        vcf_qual = '.'
        vcf_filter = 'PASS'
        vcf_info = '='.join(['RU',repeatunit]) + ';' + '='.join(['RL',str(ref_RL)])
        vcf_format = 'GT'
        vcf_sample = vcf_gt
        vcf_record = '\t'.join([str(x) for x in [chrom, vcf_start, vcf_id,
                                ref_sequence, vcf_alt, vcf_qual, vcf_filter,
                                vcf_info, vcf_format, vcf_sample] ])
        vcf_truth.write(vcf_record + '\n')

if __name__ == '__main__':
    main()

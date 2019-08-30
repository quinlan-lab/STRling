#!/usr/bin/env python

from argparse import (ArgumentParser, FileType)
import random
import sys

__author__ = 'Harriet Dashnow'
__credits__ = ['Harriet Dashnow']
__license__ = 'MIT'
__version__ = '0.1.0'
__email__ = 'h.dashnow@gmail.com'

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = ArgumentParser(description=('Generate STR genotypes where one allele is fixed and the other varies randomly in the range given.'
        'Example tab-delimited output: chr1  1  100  CAG_-2/1'))
    parser.add_argument(
        '--locus', type=str, required=True,
        help='Locus and repeat unit in the format: "chr  start  stop  repeatunit" (must be in quotes and whitespace-delimited).')
    parser.add_argument(
        '--num', type=int, default=100,
        help='Number of alleles to simulate. (default: %(default)s)')
    parser.add_argument(
        '--min', type=int, default=1,
        help='Minimum allele to simulate, in repeat units relative to the reference. (default: %(default)s)')
    parser.add_argument(
        '--max', type=int, default=400,
        help='Maximum allele to simulate, in repeat units relative to the reference. (default: %(default)s)')
    parser.add_argument(
        '--fixed', type=int, default=0,
        help='Size of fixed allele, in repeat units relative to the reference. (default: %(default)s)')
    parser.add_argument(
        '--out', type=str, default='',
        help='Prefix for output bed files.')
    parser.add_argument(
        '--seed', required=False,
        help='Random seed (can be any hashable input).')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    prefix = args.out
    assert len(args.locus.split()) == 4

    locus_string = '\t'.join(args.locus.split()) + '_{}/'.format(args.fixed)
    chrom = args.locus.split()[0]
    start = args.locus.split()[1]
    RU = args.locus.split()[3]
    allele1 = args.fixed
    if args.seed:
        random.seed(args.seed)

    # Write bed files
    for i, random_allele in enumerate(random.sample(range(args.min, args.max), args.num)):
        outfile = "{}{}.{}-{}_{}_{}_{}.bed".format(prefix, i, chrom, start, RU, allele1, random_allele)
        with open(outfile, "w") as bed_out:
            bed_out.write(locus_string + str(random_allele) + '\n')

if __name__ == '__main__':
    main()

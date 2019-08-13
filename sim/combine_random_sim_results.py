import glob
import sys
import os
import pandas as pd
from argparse import (ArgumentParser, FileType)

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = ArgumentParser(description=('Gather together simulation parameters and str results in a single file.'))
    parser.add_argument(
        '--bed_dir', type=str, default='.',
        help='Directory containing simulation design in the form of multiple bed files e.g. chr1  1  100  CAG_-2/1 (default: %(default)s)')
    parser.add_argument(
        '--str_dir', type=str, default='.',
        help='Directory containing str results. (default: %(default)s)')
    parser.add_argument(
        '--out', type=str, default='',
        help='Prefix for output csv files. (default: %(default)s)')
    return parser.parse_args()

def get_sim(f):
    return(os.path.basename(f).split('.')[0])

def get_sim_str(f):
    return(os.path.basename(f).split('.')[0].split('_')[0])

def parse_alleles(s):
    repeatunit = s.split('_')[0]
    sizes = s.split('_')[1].split('/')
    return [repeatunit] + sizes

def parse_bed(f):
    with open(f) as bed:
        firstline = bed.readlines()[0].split()
        alleles = parse_alleles(firstline[-1])
        return(alleles)

def parse_bounds(all_files):
    header = 'chrom, left, right, mean-position-of-non-split reads, number-of-left-splits, number-of-right-splits, n-total-str-reads, repeat-unit'
    all_df_bounds = []
    for f in all_files:
        df = pd.read_csv(f, sep='\t', names=header.split(','))
        df['sim'] = get_sim_str(f)
        all_df_bounds.append(df)
    all_df = pd.concat(all_df_bounds)
    return(all_df)

def main():
    # Parse command line arguments
    args = parse_args()

    bedfiles = glob.glob(args.bed_dir+"/*.bed")
    if len(bedfiles) == 0:
        sys.exit('ERROR: No bedfiles found in the given directory: ' + args.bed_dir)
    bed_df = pd.DataFrame(columns=['repeatunit', 'allele1', 'allele2', 'sim'])
    for f in bedfiles:
        sim = get_sim(f)
        bed_df.loc[int(sim)] = parse_bed(f)+[sim]
    print(bed_df)
    bed_df.to_csv(args.out + '-sims.csv', index=False)
    
    boundsfiles = glob.glob(args.str_dir+"/*-bounds.txt")
    if len(boundsfiles) == 0:
        sys.exit('ERROR: No -bounds.txt files found in the given directory: ' + args.str_dir)
    parse_bounds(boundsfiles).to_csv(args.out + '-bounds.csv', index=False)

    readsfiles = glob.glob(args.str_dir+"/*-reads.txt")
    if len(readsfiles) == 0:
        sys.exit('ERROR: No -reads.txt files found in the given directory: ' + args.str_dir)
    all_df_reads = []
    for f in readsfiles:
        df_reads = pd.read_csv(f, sep='\t')
        df_reads['sim'] = get_sim_str(f)
        all_df_reads.append(df_reads)
    all_df = pd.concat(all_df_reads)
    all_df = all_df.merge(bed_df, how='left', on='sim')
    all_df.to_csv(args.out + '-results.csv', index=False)

if __name__ == "__main__":
    main()

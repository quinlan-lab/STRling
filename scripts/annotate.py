# Parse TRF ngs dat file into tab delimited file
# Assumes TRF was run on aligned contigs.
# Each @ is a contig labeled with it's mapping position,
# each line is an STR

import argparse
import sys
import os
from strtools import normalise_str
import pandas as pd
import pyranges as pr

def parse_args():
    # top-level parser
    parser = argparse.ArgumentParser()

    parser.add_argument('--strling', type=str,
                        help='STRling outlier tsv file for one sample')
    parser.add_argument('--known', type=str,
                        help='known STRs in bed format with a 4th column containing the repeat unit/motif')
    parser.add_argument('--centromeres', type=str, help='Optional bed file of centromeres to annotate')
    parser.add_argument('--telomeres', type=str, help='Optional bed file of telomeres to anotate')
    parser.add_argument('--LCRs', type=str, help='Optional bed file of Low copy repeats to anotate')
    parser.add_argument('--slop', type=int, default=200,
                        help='Consider loci the same if they are no more than this distance apart and have the same repeat unit')
    parser.add_argument('--out', type=str, default = '',
                        help='output file of variants in annotated bed format')

    return parser.parse_args()

def parse_str_bed(filename):
    """Parse bed file with columns chr,start,end,repeatunit
        and optional header row starting with #"""
    try:
        df = pd.read_csv(filename, delim_whitespace = True)
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    df.columns = ['Chromosome', 'Start', 'End', 'repeatunit']

    for repeatunit in df['repeatunit']:
        for base in repeatunit:
            if base not in ['A', 'T', 'C', 'G']:
                sys.exit('ERROR: Non-DNA found in the third column of {}: {}\n'.format(filename, repeatunit))

    return(df)

def match_closest(strling_df, known_df, this_repeatunit, slop):
    """Filter to a specific repeat unit then annotate with the closest locus"""
    this_strling_df = strling_df.loc[strling_df['repeatunit_norm'] == this_repeatunit].copy()

    if this_strling_df.empty:
        return this_strling_df

    this_known_df = known_df[known_df['repeatunit_norm'] == this_repeatunit]
    if this_strling_df.empty:
        return this_strling_df #XXX May need to add extra columns?

    strling_pr = pr.PyRanges(this_strling_df)
    known_pr = pr.PyRanges(this_known_df)

    # Annotate with the closest locus
    nearest_pr = strling_pr.nearest(known_pr)
    nearest_df = nearest_pr.df

    if nearest_df.empty:
        return this_strling_df
    #print("before")
    #print(this_strling_df)

    # Remove reference STRs more than slop bp away
    nearest_columns = ['Start_b', 'End_b', 'repeatunit_norm_b', 'Distance']
    #nearest_df.loc[nearest_df.Distance > slop, nearest_columns] = None
    nearest_df['annotated'] = None
    nearest_df.loc[nearest_df.Distance <= slop, 'annotated'] = 'reference'
    nearest_df.loc[nearest_df.Distance > slop, 'annotated'] = 'novel'
    nearest_df['annotated'].fillna('novel')

    # Create a unique index
    this_strling_df['locus'] = this_strling_df['Chromosome'] + '-' + this_strling_df['Start'].astype(str
        ) + '-' + this_strling_df['End'].astype(str) + '-' + this_strling_df['repeatunit']
    this_strling_df.set_index('locus', inplace = True)
    nearest_df['locus'] = nearest_df['Chromosome'].astype(str) + '-' + nearest_df['Start'
        ].astype(str) + '-' + nearest_df['End'].astype(str) + '-' + nearest_df['repeatunit']
    nearest_df.set_index('locus', inplace = True)

    #print('nearest')
    #print(nearest_df)

    nearest_df = nearest_df.filter(['annotated', 'Distance'])
    this_strling_df = this_strling_df.merge(nearest_df, how = 'left', left_index = True,
                                    right_index = True)
    #print('after')
    #print(this_strling_df)
    #print('\n')
    return this_strling_df

def match_variants(strling_df, known_df, slop):
    """Match known STRs with strling variants
        - Within X bp of slop
        - Same repeat unit"""

    strling_df['repeatunit_norm'] = strling_df['repeatunit'].apply(normalise_str)
    known_df['repeatunit_norm'] = known_df['repeatunit'].apply(normalise_str)

    # Break down by repeat unit before comparing
    all_closest_df = pd.DataFrame()
    for this_ru in set(strling_df['repeatunit_norm']):
        all_closest_df = all_closest_df.append(match_closest(strling_df, known_df, this_ru, slop))
    
    return all_closest_df

def annotate_feature_cov(strling_df, bed, feature):
    """Annotate strling calls with if they overlap a feature"""
    strling_pr = pr.PyRanges(strling_df)
    bed_pr = pr.readers.read_bed(bed)

    cov_pr = strling_pr.coverage(bed_pr)
    cov_df = cov_pr.df

    # Create a unique index
    cov_df['locus'] = cov_df['Chromosome'].astype(str) + '-' + cov_df['Start'
        ].astype(str) + '-' + cov_df['End'].astype(str) + '-' + cov_df['repeatunit']
    cov_df.set_index('locus', inplace = True)

    cov_df = cov_df.filter(['NumberOverlaps'])
    cov_df.rename(columns={'NumberOverlaps': f'{feature}Cov'}, inplace = True)
    strling_df = strling_df.merge(cov_df, how = 'left', left_index = True,
                                    right_index = True)
    strling_df[f'{feature}Cov'] = strling_df[f'{feature}Cov'].fillna(0)

    return(strling_df)

def main():
    args = parse_args()

    # Parse inputs as pandas data frame (df) or pyranges (pr) objects
    strling_df = pd.read_csv(args.strling, sep='\t')
    strling_df.dropna(subset = ['chrom'], inplace = True) # Remove rows where chromosome is missing
    strling_df = strling_df.rename(columns={'chrom': 'Chromosome',
                    'left': 'Start', 'right': 'End'})
    # +1 end of any regions of length 0 (Start == End) to avoid pyranges error
    strling_df.loc[strling_df['Start'] == strling_df['End'], 'End'] += 1

    known_df = parse_str_bed(args.known)

    # Annotate strling calls with known STRs
    strling_df = match_variants(strling_df, known_df, args.slop)

    # Annotate with additional stuff if requested
    if args.centromeres:
        strling_df = annotate_feature_cov(strling_df, args.centromeres, 'centromere')
    if args.telomeres:
        strling_df = annotate_feature_cov(strling_df, args.telomeres, 'telomere')
    if args.LCRs:
        strling_df = annotate_feature_cov(strling_df, args.LCRs, 'LCR')

    strling_df.to_csv(args.out, sep='\t', index=False)

if __name__ == '__main__':
    main()

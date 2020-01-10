#!/usr/bin/env python
"""Read STRling output and look for individuals that are outliers at STR loci
"""

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)
    import time
    import datetime
    import argparse
    import sys
    import glob
    import os
    #import re
    import numpy as np
    import statsmodels.api as sm
    from scipy.stats import norm
    from statsmodels.sandbox.stats.multicomp import multipletests
    import pandas as pd
    import pyranges as pr

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Read STRling output and look for individuals that are outliers at STR loci')
    parser.add_argument(
        '--genotypes', type=str, nargs='+', required = True,
        help='-genotype.txt files for all samples produced by STRling.')
    parser.add_argument(
        '--unplaced', type=str, nargs='+', required = True,
        help='-unplaced.txt files for all samples produced by STRling. Contains the number of unassigned STR reads for each repeat unit.')
    parser.add_argument(
        '--out', type=str, default = '',
        help='Prefix for all output files (suffix will be STRs.tsv) (default: %(default)s)')
    parser.add_argument(
        '--control', type=str, default='',
        help='Input file for median and standard deviation estimates at each locus from a set of control samples. This file can be produced by this script using the emit option. If this option is not set, all samples in the current batch will be used as controls by default.')
    parser.add_argument(
        '--emit', type=str, default='',
        help='Output file for median and standard deviation estimates at each locus (tsv).')
    parser.add_argument(
        '--slop', type=int, default=50,
        help='Merge loci that are within this many bp of each other and have the same repeat unit.')
    parser.add_argument(
        '--min_clips', type=int, default=0,
        help='In the individual sample files, only report loci with at least many soft-cliped reads in that sample.')
    parser.add_argument(
        '--min_size', type=int, default=0,
        help='In the individual sample files, only report loci with at least this allele2_est size in that sample.')

    return parser.parse_args()

def convert_time(s):
    """Convert time in seconds to time in hours:minutes:seconds"""
    return str(datetime.timedelta(seconds = int(s)))

def get_sample(fullpath):
    """Get the sample ID from the filename"""
    basename = os.path.basename(fullpath)
    return(basename.rsplit('-', maxsplit = 1)[0])

def parse_unplaced(filename):
    """Parse unplaced STR read counts"""
    sample_id = get_sample(filename)
    try:
        unplaced_counts = pd.read_csv(filename, delim_whitespace = True, header = None,
                            names = ['repeatunit', 'unplaced_count'])
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    unplaced_counts['sample'] = sample_id
    unplaced_counts = unplaced_counts[['sample', 'repeatunit', 'unplaced_count']]
    return(unplaced_counts)

def parse_genotypes(filename, min_clips = 5):
    """Parse -genotype.txt file produced by STRling"""
    sample_id = get_sample(filename)
    try:
        genotype_data = pd.read_csv(filename, delim_whitespace = True, header = 0)
        genotype_data.rename(columns={'#chrom': 'Chromosome', 'left': 'Start', 'right': 'End'}, inplace = True)
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    if genotype_data.shape[0] == 0: # Check for file with only header
        sys.exit('ERROR: file {0} contained 0 loci.\n'.format(filename))

    # Report number of samples
    sys.stderr.write('Sample: {} Loci: {}\n'.format(sample_id, genotype_data.shape[0]))

    genotype_data['sample'] = sample_id
    return(genotype_data)

def parse_controls(control_file):
    """Parse control file with columns locus, median and standard deviation"""

    control_estimates = pd.read_csv(control_file, index_col=0, delim_whitespace = True, header = None)

    # Allow for old style column headings, but change to mu and sd.
    if control_estimates.columns[0] in ['mu', 'median'] and control_estimates.columns[1] in ['sd', 'SD']:
        colnames = list(control_estimates.columns)
        colnames[0:2] = ['mu', 'sd']
        control_estimates.columns = colnames
    else:
        raise ValueError(''.join(["The column names in the control file ",
        "don't look right, expecting columns named median, SD ",
        "or mu, sd. Column names are ", str(list(control_estimates.columns)),
        ". Check the file: ", control_file]))
    return(control_estimates)

def z_score(x, df):
    """Calculate a z score for each x value, using estimates from a pandas data
    frame with the columns 'mu' and 'sd' and index coressponding to the x values"""
    z = (x.transpose() - df['mu'])/df['sd']
    return z.transpose()

def p_adj_bh(x):
    '''Adjust p values using Benjamini/Hochberg method'''
    # Mask out nan values as they cause the multiptests algorithm to return all nan
    mask = np.isfinite(x)
    pval_corrected = x
    pval_corrected[mask] = multipletests(x[mask], method='fdr_bh', returnsorted = False)[1]
    return pval_corrected

def main():

    start_time = time.time()

    # Parse command line arguments
    args = parse_args()

    base_filename = args.out
    emit_file = args.emit
    control_file = args.control
    genotype_files = args.genotypes
    unplaced_files = args.unplaced
    slop = args.slop
    min_clips = args.min_clips
    min_size = args.min_size
    results_suffix = 'STRs.tsv'

    # Check files exist for all samples
    genotype_ids = set([get_sample(f) for f in genotype_files])
    unplaced_ids = set([get_sample(f) for f in unplaced_files])
    if genotype_ids == unplaced_ids:
        all_samples = genotype_ids
    else:
        all_samples = genotype_ids | unplaced_ids
        missing_samples = (all_samples - genotype_ids) | (all_samples - unplaced_ids)
        sys.exit("ERROR: One or more files are missing for sample(s): " + ' '.join(missing_samples))
    
    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Reading input files for {0} samples\n'.format(len(all_samples)))

    if len(all_samples) < 2 and control_file == '':
        sys.stderr.write('WARNING: Only 1 sample and no control file provided, so outlier scores and p-values will not be generated.')

    # Parse unplaced data
    unplaced_data = pd.concat( (parse_unplaced(f) for f in unplaced_files), ignore_index = True)

    unplaced_wide = unplaced_data
    # Fill zeros in unplaced counts
    unplaced_wide = unplaced_data.pivot(index='repeatunit', columns='sample',
                    values='unplaced_count').fillna(0)
    unplaced_wide['repeatunit'] = unplaced_wide.index

    sample_cols = list(set(unplaced_data['sample']))
    unplaced_long = pd.melt(unplaced_wide, id_vars = 'repeatunit',
                            value_vars = sample_cols, value_name = 'unplaced_count',
                            var_name = 'sample')

    # Write unplaced read counts
    unplaced_filename = base_filename + 'unplaced.tsv'
    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Writing unplaced counts for all samples to {}\n'.format(unplaced_filename))
    unplaced_long.to_csv(unplaced_filename, sep= '\t', index = False, na_rep='NaN')

    # Parse genotype data
    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Reading genotypes\n')
    genotype_data = pd.concat( (parse_genotypes(f) for f in genotype_files), ignore_index = True)

    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Processing {} samples...\n'.format(len(all_samples)))

    #XXX fix column names in parse_genptypes()
    genotype_data['locus'] = genotype_data['Chromosome'] + '-' + genotype_data['Start'].astype(str) + '-' + genotype_data['End'].astype(str) + '-' + genotype_data['repeatunit']
    genotype_data = genotype_data.rename(columns={'Chromosome': 'chrom', 'Start': 'left', 'End': 'right'})
    genotype_data = genotype_data.reset_index(drop=True)

#XXX some testing here - remove
#    genotype_data.reset_index(drop=True, inplace=True)
#    genotype_data.to_csv("test.csv")
#    genotype_data.loc[genotype_data.duplicated(subset = ['locus', 'sample'])].to_csv("test-duplicates.tsv", sep = '\t')

    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Calculating median sample depths\n')
    # Calculate median depth per sample
    sample_depths = genotype_data[['sample','depth']].groupby('sample').median(skipna=True)
    sample_depths['sample'] = sample_depths.index
    sample_depths.to_csv(base_filename + 'depths.tsv', sep= '\t', index = False, na_rep='NaN')

    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Fill zeros\n')
    # Fill NA in genotype
    # Not appropriate to fill with 0 because NA are likely regions of very high coverage 
    sum_str_wide = genotype_data.pivot(index='locus', columns='sample',
                    values='sum_str_counts')

    sample_cols = list(set(genotype_data['sample']))

    # Remove rows that are all 0 or all NA
    #sum_str_wide[(sum_str_wide.sum(axis = ) == 0) | (x[:,1] == 0.)]
    mask = np.all(np.isnan(sum_str_wide) | np.equal(sum_str_wide, 0), axis=1)
    sum_str_wide = sum_str_wide[~mask]
    sum_str_wide['locus'] = sum_str_wide.index

    sum_str_long = pd.melt(sum_str_wide, id_vars = 'locus',
                            value_vars = sample_cols, value_name = 'sum_str_counts',
                            var_name = 'sample')

    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Add locus info back in and fill additional zeros\n')
    # Add locus info back in 
    genotype_data = pd.merge(sum_str_long, genotype_data, how='left')
    # Fill zeros in additional columns
    genotype_data[['left', 'right']] = genotype_data[['left',
                                                    'right']].fillna(0)

    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Fill median sample depths\n')
    # Use the median sample depth as a proxy for local depth for loci with zero or NA depth
    genotype_data['depth'] = genotype_data['depth'].replace({0: np.nan})
    genotype_data['depth'] = genotype_data[['sample','depth']].groupby('sample'
                            ).transform(lambda x: x.fillna(x.median(skipna=True)))


    # Normalise STR coverage by median coverage
    factor = 1 # This was 100 in STRetch

    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Log normalize\n')
    genotype_data['sum_str_log'] = np.log2(factor * (genotype_data['sum_str_counts'] + 1) / genotype_data['depth'])
   
    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Pivot data for z score calculations\n')
    # For each locus, calculate if that sample is an outlier relative to others

    sample_depths = genotype_data[['sample', 'depth']].groupby('sample').median()
    # Calculate values for if there were zero reads at a locus in all samples
    null_locus_counts = np.log2(factor * (0 + 1) / sample_depths)
    # Add a null locus that has 0 reads for all individuals (so just uses coverage)
    null_locus_counts_est = pd.Series([
        np.median(null_locus_counts['depth']), np.std(null_locus_counts['depth'])
        ], index = ['mu', 'sd'])
    #XXX do something if sd == 0?

    # Calculate a z scores using median and SD estimates from the current set
    # of samples

    # Use Huber's M-estimator to calculate median and SD across all samples
    # for each locus
    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Calculate mu and sd estimates\n')
    sum_str_wide.drop(['locus'], axis = 1, inplace = True)
    locus_estimates = pd.DataFrame(list(zip(
        np.median(sum_str_wide, axis = 1), np.std(sum_str_wide, axis = 1)
        )), columns = ['mu', 'sd'])
    locus_estimates.index = sum_str_wide.index

    # Where sd is NA, replace with the minimum non-zero sd from all loci
    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Clean up zeros in sd estimates\n')
    min_sd = np.min(locus_estimates['sd'][locus_estimates['sd'] > 0])
    locus_estimates['sd'].fillna(min_sd, inplace=True)
    # if sd is 0, replace with min_sd #XXX is this sensible?
    if null_locus_counts_est['sd'] == 0:
        null_locus_counts_est['sd'] = min_sd

    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Export and/or import control files, if any\n')
    # Save median and SD of all loci to file if requested (for use as a
    # control set for future data sets)
    if emit_file != '':

        locus_estimates.loc['null_locus_counts'] = null_locus_counts_est

        n = len(sum_str_wide.columns)
        locus_estimates['n'] = n

        locus_estimates.to_csv(emit_file, sep= '\t')

    # Calculate z scores using median and SD estimates per locus from a
    # provided control set
    if control_file != '':
        # Parse control file
        control_estimates = parse_controls(control_file)
        # Get a list of all loci in the control file but not the sample data
        control_loci_df = control_estimates.iloc[control_estimates.index != 'null_locus_counts']
        control_loci = [x for x in control_loci_df.index if x not in sum_str_wide.index]

        # Extract and order just those control estimates appearing in the current data
        mu_sd_estimates = control_estimates.reindex(sum_str_wide.index)
        # Fill NaNs with null_locus_counts values
        mu_sd_estimates.fillna(control_estimates.loc['null_locus_counts'],
                                inplace=True)
    else:
        # Extract and order estimates to match the current data
        mu_sd_estimates = locus_estimates.reindex(sum_str_wide.index)

    # calculate z scores
    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Calculating z scores\n')
    z = z_score(sum_str_wide, mu_sd_estimates)

    # If a control file is given, effectively add zeros counts at all loci in 
    # controls but not in the samples. 
    # These extra rows will dissapear due to a later merge
    if control_file != '': 
        # Create a sum_str_wide as if all loci have zero counts
        null_sum_str_wide = pd.DataFrame(columns = sample_names, index = control_loci)
        null_sum_str_wide.fillna(null_locus_counts, inplace = True)
        # Caculate z scores
        null_z = z_score(null_sum_str_wide, 
                            control_estimates.reindex(null_sum_str_wide.index))
        loci_with_counts = z.index
        z = z.append(null_z)

    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Calculating and adjusting p values\n')
    if z.shape[0] == 1:
        ids = z.columns # save index order as data gets sorted
        # Calculate p values based on z scores (one sided)
        z_list = list(z.iloc[0])
        pvals = norm.sf(z_list) # no need to adjust p values if one locus

        # Merge pvals and z scores back into genotype_data
        p_z_df = pd.DataFrame({'sample': ids, 'p_adj': pvals, 'outlier': z_list})
        genotype_data = pd.merge(genotype_data, p_z_df)

    elif z.shape[0] > 1:
        # Calculate p values based on z scores (one sided)
        pvals = z.apply(lambda z_row: [norm.sf(x) for x in z_row], axis=1, 
                    result_type='broadcast') # apply to each row
        if pvals.isnull().values.all(): # Don't bother adjusting p values if all are null
            adj_pvals = pvals
        else:
            # Adjust p values using Benjamini/Hochberg method
            adj_pvals = pvals.apply(p_adj_bh, axis=0) # apply to each column
        
        # Merge pvals and z scores back into genotypes
        adj_pvals['locus'] = adj_pvals.index
        adj_pvals_long = pd.melt(adj_pvals, id_vars = 'locus',
                                value_vars = sample_cols, value_name = 'p_adj',
                                var_name = 'sample')
        genotype_data = pd.merge(genotype_data, adj_pvals_long)

        pvals['locus'] = pvals.index
        pvals_long = pd.melt(pvals, id_vars = 'locus',
                                value_vars = sample_cols, value_name = 'p',
                                var_name = 'sample')
        genotype_data = pd.merge(genotype_data, pvals_long)

        z['locus'] = z.index #important to do this only after p values calculated
        z_long = pd.melt(z, id_vars = 'locus',
                        value_vars = sample_cols, value_name = 'outlier', var_name = 'sample')
        genotype_data = pd.merge(genotype_data, z_long)

    elif z.shape[0] == 0:
        raise ValueError('z score table is empty')

    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Write output files\n')
    # Specify output data columns
    write_data = genotype_data[['chrom', 'left', 'right', 'locus',
                                    'sample', 'repeatunit',
                                    'allele1_est', 'allele2_est',
                                    'spanning_reads', 'spanning_pairs',
                                    'left_clips', 'right_clips', 'unplaced_pairs',
                                    'sum_str_counts', 'sum_str_log', 'depth',
                                    'outlier', 'p', 'p_adj',
                                    ]]

    #sort by outlier score then estimated size (bpInsertion), both descending
    write_data = write_data.sort_values(['p_adj', 'allele2_est'], ascending=[True, False])
    # Convert outlier and p_adj to numeric type and do some rounding/formatting
    write_data['outlier'] = pd.to_numeric(write_data['outlier'])
    write_data['p_adj'] = [ format(x, '.2g') for x in pd.to_numeric(write_data['p_adj']) ]
    write_data = write_data.round({'outlier': 1, 'sum_str_log': 1,
                                    'sum_str_log': 1})
    int_cols = ['left', 'right', 'sum_str_counts', 'spanning_reads', 'spanning_pairs',
                'left_clips', 'right_clips', 'unplaced_pairs']
    write_data[int_cols] = write_data[int_cols].astype('Int64')

    # Write individual files for each sample, remove rows where locuscoverage == 0
    samples = set(write_data['sample'])
    for sample in samples:
        sample_filename = base_filename + sample + '.' + results_suffix
        sample_df = write_data.loc[write_data['sample'] == sample]
        # filter to loci with a minimum number of clipped reads in that sample XXX change?
        sample_df = sample_df[sample_df['allele2_est'] >= min_size]
        sample_df = sample_df[sample_df['left_clips'] + sample_df['right_clips'] >= min_clips]
        sample_df.to_csv(sample_filename, sep= '\t', index = False, na_rep='NaN')

    # Write all samples to a single file
    all_filename = base_filename + results_suffix
    write_data.to_csv(all_filename, sep= '\t', index = False, na_rep='NaN')
    sys.stderr.write(f'Elapsed time: {convert_time(time.time() - start_time)} ')
    sys.stderr.write('Finished!\n')

if __name__ == '__main__':
    main()

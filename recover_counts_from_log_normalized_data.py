import h5py
import numpy as np
import pandas as pd
import sys
from optparse import OptionParser
import json

# Small number for checking integer values in presence of floating point error
EPSILON = 1e-10

# Maximum default size-factor search range
MAX_RANGE = 1e5

def main():
    usage = "%prog [options] input_file" 
    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--data_type", help="Data type (required). Must be one of: 'TSV', 'CSV'")
    parser.add_option("-l", "--log_base", help="Base of logarithm used for normalization. If not provided, the natural logarithm is assumed.")
    parser.add_option("-m", "--mult_factor", help="Multiplicative factor used for normalized (e.g., for TPM, this is one million)")
    parser.add_option("-a", "--max_size_factor", help="Maximum size-factor search value for each cell/sample. Default: 100k (for 10x scRNA-seq data)")
    parser.add_option("-t", "--transpose", action='store_true', help="Transposed matrix. If True, the input matrix is gene-by-cell. If False, the input matrix is cell-by-gene")
    parser.add_option("-v", "--verbose", action='store_true', help="Verbose. If True, output logging information.")
    parser.add_option("-o", "--output_file", help="Output file")
    (options, args) = parser.parse_args()

    # Input file
    in_f = args[0]
    
    # Read parameters
    if options.log_base:
        log_base = float(options.log_base)
    else:
        log_base = None
    
    if options.mult_factor:
        mult_factor = float(options.mult_factor)
    else:
        raise Exception("Please provide a multiplicative factor using the '-m' option.")

    out_f = options.output_file

    if options.data_type == 'TSV':
        sep = '\t'
    elif options.data_type == 'CSV':
        sep = ','
    else:
        raise Exception("Please provide a data type (TSV or CSV) using the '-d' option.")

    if options.max_size_factor:
        max_range = int(options.max_size_factor)
    else:
        max_range = MAX_RANGE

    verbose = options.verbose

    df_expr = pd.read_csv(in_f, sep=sep, index_col=0)

    # Get cell-by-gene expression matrix
    X = np.array(df_expr)
    if options.transpose:
        X = X.T

    # Recover counts
    counts, size_factors = recover_counts(
        X, 
        mult_factor,
        max_range, 
        log_base=log_base, 
        verbose=verbose
    )
    if options.transpose:
        counts = counts.T

    print(size_factors)
    print(np.median(size_factors))
    
    # Create output dataframe
    if verbose:
        print(f"Writing output to {out_f}...")
    df_out = pd.DataFrame(
        data=counts,
        columns=df_expr.columns,
        index=df_expr.index
    )

    # Write output
    df_out.to_csv(out_f, sep=sep)


def recover_counts(X, mult_value, max_range, log_base=None, verbose=True):
    """
    Given log-normalized gene expression data, recover the raw read/UMI counts by
    inferring the unknown size factors.

    Parameters
    ----------
    X: 
        The log-normalized expression data. This data is assumed to be normalized 
        via X := log(X/S * mult_value + 1)
    max_range:
        Maximum size-factor search range to use in binary search.
    mult_value:
        The multiplicative value used in the normalization. For example, for TPM
        this value is one millsion. For logT10K, this value is ten thousand.
    log_base:
        The base of the logarithm

    Returns
    -------
    counts:
        The inferred counts matrix
    size_factors:
        The array of inferred size-factors (i.e., total counts)
    """
    # Recover size-factor normalized ata
    if log_base is None:
        X = (np.exp(X) - 1) / mult_value
    else:
        X = (np.pow(X, log_base)-1) / mult_value

    # Infer the size factors
    size_factors = []
    for x_i, x in enumerate(X):
        if verbose and x_i % 100 == 0:
            print(f"Found size factor in {x_i} cells...")
        s = binary_search(x, max_r=max_range)
        size_factors.append(s)
    counts = X.T * size_factors
    counts = (counts.T).astype(int)
    return counts, size_factors


def binary_search(vec, min_r=0, max_r=100000):
    """
    For a given gene expression vector corresponding to a single cell
    or gene expression profile, infer the unknown size factor using
    binary search.
    """
    # Find the smallest non-zero expression value. We assume
    # that this value corresponds to a count of one.
    min_nonzero = np.min([x for x in vec if x != 0])

    # Initialize the search bounds
    min_bound = min_r
    max_bound = max_r

    # Run binary search
    while True:
        curr_s = int((max_bound - min_bound) / 2) + min_bound
        cand_count = min_nonzero * curr_s

        if np.abs(cand_count - 1) < EPSILON:
            return curr_s
        elif cand_count > 1: # The size factor is too big
            max_bound = curr_s 
        elif cand_count < 1: # The size factor is too small
            min_bound = curr_s

        # Sometimes the floating point percision is higher than our check. Instead of relaxing 
        # our tolerance, we simply choose the size-factor with the smaller difference between
        # the putative one-count and the true putative one-count.
        if max_bound - min_bound == 1:
            cand_count_max = min_nonzero * max_bound
            cand_count_min = min_nonzero * min_bound
            diff_max = np.abs(cand_count_max - 1)
            diff_min = np.abs(cand_count_min - 1)
            if diff_max < diff_min:
                return max_bound
            else:
                return min_bound

if __name__ == "__main__":
    main()

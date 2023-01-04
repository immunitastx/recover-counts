# Recover counts from log-normalized scRNA-seq data

This script recovers raw counts from log-normalized gene expression data. Specifically, this function assumes that
each cell/sample was normalized via: X = log(C/S * M + 1) where C are the counts, S are the total counts in that 
cell/sample, and M is a multiplicative factor (e.g., for TPM data, this value would be one million). Given X, this 
program recovers C. 

This algorithm works as follows: first, we assume that the smallest non-zero count in each cell/sample is 1. Thus, if we 
find this value, and we know the correct value for S, then the inverse of the normalization function should equal one. 
That is, 1 = (exp(X) - 1) * (S/M). Unfortunately, we don't know S. This algorithm  performs a binary search over possible 
values for S, the total counts, until it finds the value for which we can recover a count of one using the inverse of the 
log-normalization function.

```
Usage: recover_counts_from_log_normalized_data.py [options] input_file

Options:
  -h, --help            show this help message and exit
  -d DATA_TYPE, --data_type=DATA_TYPE
                        Data type (required). Must be one of: 'TSV', 'CSV'
  -l LOG_BASE, --log_base=LOG_BASE
                        Base of logarithm used for normalization. If not
                        provided, the natural logarithm is assumed.
  -m MULT_FACTOR, --mult_factor=MULT_FACTOR
                        Multiplicative factor used for normalized (e.g., for
                        TPM, this is one million)
  -a MAX_COUNT, --max_count=MAX_COUNT
                        Maximum total UMI count/read depth for each
                        cell/sample. Default: 100k (for 10x scRNA-seq data)
  -t, --transpose       Transposed matrix. If True, the input matrix is gene-
                        by-cell. If False, the input matrix is cell-by-gene
  -v, --verbose         Verbose. If True, output logging information.
  -o OUTPUT_FILE, --output_file=OUTPUT_FILE
                        Output file
```

## Unit testing

To test the code, download a test scRNA-seq dataset from GEO (from [Laughney et al.](https://doi.org/10.1038/s41591-019-0750-6)) via the following command:

`curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3516nnn/GSM3516673/suppl/GSM3516673_MSK_LX682_NORMAL_dense.csv.gz`

Then, run the unit test via:

`python -m unittest test.py`

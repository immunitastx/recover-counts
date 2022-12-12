# Recover counts from log-normalized scRNA-seq data

This script recovers raw counts from log-normalized gene expressiond data. Specifically, this function assumes that
each cell/sample was normalized via: X = log(C/S * M + 1) where C are the counts, S are the total counts in that 
cell/sample, and M is a multiplicative factor (e.g., for TPM data, this value would be one million). 

Note this function assumes that the smallest non-zero count in each cell/sample is 1 and performs a binary search over possible
size-factors until it finds the size factor for which we can recover a count of one using the inverse of the log-normalization
function.

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
  -a MAX_SIZE_FACTOR, --max_size_factor=MAX_SIZE_FACTOR
                        Maximum size-factor search value for each cell/sample.
                        Default: 100k (for 10x scRNA-seq data)
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

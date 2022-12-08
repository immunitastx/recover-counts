# Recover counts from log-normalized scRNA-seq data

This script recovers raw counts from log-normalized gene expressiond data. Specifically, this function assumes that
each cell/sample was normalized via: X = log(C/S * M + 1) where C are the counts, S are the total counts in that 
cell/sample, and M is a multiplicative factor (e.g., for TPM data, this value would be one million). 


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
  -t, --transpose       Transposed matrix. If True, the input matrix is gene-
                        by-cell. If False, the input matrix is cell-by-gene
  -v, --verbose         Verbose. If True, output logging information.
  -o OUTPUT_FILE, --output_file=OUTPUT_FILE
                        Output file
```

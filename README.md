# PLSDA-multiple
# Language: R
# Input: TXT (parameters)
# Output: prefix
# Tested with: PluMA 1.0, R 3.2.5

PluMA plugin that runs Partial Least Squares Discriminant Analysis (PLS-DA, Stahle and Wold 1987) 
to extract features from a set of viral samples, for the use of downstream analysis.

The input file will be a textfile with each line consisting of tab-deliminted keyword-value
pairs.  The only parameters for this plugin are: "training" and "clinical" data.

The outputfile is a prefix, all chosen features will be output in the following format:
outputfile_(time).csv

Note the CSV files in the example/ directory are not publically available.
A future goal is to make a synthetic data set.  In the meantime however, one may
use it on their own input data.

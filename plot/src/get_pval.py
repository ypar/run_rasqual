#!/usr/bin/env python

###
# YoSon
# 08/20/2015
# @ run_rasqual/plot/get_pval.py
###

import csv, sys, os
import numpy as np
from scipy.stats import chisqprob


if __name__ == '__main__':
    '''takes rasqual output file and append a column containing chisq pvalues'''


infile = str(sys.argv[1])

outfile = str(sys.argv[2])



with open(infile, 'r') as f, open(outfile, 'w') as ofile:
        
    heads = ['feature_id', 'rs_id', 'chrom', 'pos', 'ref', 'alt', 'allele_freq', 'hwe_chisq', 'impute_qual', 'impute_qua_rsq', 'chisq', 'effect_size', 'seq_err_rate_delta', 'ref_bias_phi', 'overdisp', 'snp_id_region', 'no_feature_snps', 'no_tested_snps', 'no_iter_h0', 'no_iter_h1', 'rand_loc_ties', 'likelihood_h0', 'convergence', 'r2_fsnps', 'r2_rsnps', 'chisq_pval']
        
    cout = csv.writer(ofile, delimiter='\t')
    cout.writerow(heads)
        
        
    for line in f:
        line = line.rstrip().split()
        c2stat = line[10]
            
        pval = chisqprob(float(c2stat), 1)
            
        nl = np.append(line, pval)
        cout.writerow(nl)
            


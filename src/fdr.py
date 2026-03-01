'''
False Discovery Rate (FDR) correction for phasing p-values.

Applies Benjamini-Hochberg FDR correction to the two p-value columns
(pval and pval_b) produced by the phasing statistics stage, adding
FDR-adjusted columns to the summary file in place.

Author: Chao-Jung Wu
Date:   2024-Mar-25
'''
import statsmodels.stats.multitest as mu
import pandas as pd

def cal_fdr(infile):
  """Apply FDR correction to pval and pval_b columns and overwrite *infile*."""
  df = pd.read_csv(infile, sep='\t', low_memory=False)

  pvals = df['pval'].tolist()
  fdrs = mu.fdrcorrection(pvals)
  df['pval_fdr'], df['pval_accept'] = fdrs[1], fdrs[0]

  pvals = df['pval_b'].tolist()
  fdrs = mu.fdrcorrection(pvals)
  df['pvalb_fdr'], df['pvalb_accept'] = fdrs[1], fdrs[0]

  df.to_csv(infile, sep='\t', index=False)

def demo():
  infile = 'pvalues.txt'
  with open(infile, 'r') as fh:
    pvals = [float(x.rstrip('\n')) for x in fh.readlines()]

  fdrs = mu.fdrcorrection(pvals)
  outfile = 'pvalues_fdr_accept.txt'
  with open(outfile, 'w') as fho:
    for i in range(len(pvals)):
      data = [pvals[i], fdrs[1][i], fdrs[0][i]]
      print('\t'.join(str(j) for j in data), file=fho)

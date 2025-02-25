'''
Date: 2024-03-25
Chao-Jung Wu

Say I have 1000 contigs.
Once I obtain the 1000 pvalues for each of the 1000 contigs, 
FDR determines which pvalue is indeed significant (p < 0.05) against the pvalue population at level alpha (default 0.05).
'''
import statsmodels.stats.multitest as mu
import pandas as pd

def demo ():
  infile = 'pvalues.txt'
  with open (infile, 'r') as fh:
    pvals = [float(x.rstrip('\n')) for x in fh.readlines()]
  
  fdrs = mu.fdrcorrection(pvals)
  #for i in fdrs[0]: print(i) # TURE accept; FALSE reject
  #for i in fdrs[1]: print(i) # FDR, adjusted pvalue

  outfile = 'pvalues_fdr_accept.txt'
  fho = open (outfile, 'w')
  for i in range(len(pvals)):
    data = [pvals[i], fdrs[1][i], fdrs[0][i]]
    print('\t'.join(str(j) for j in data), file=fho)
  fho.close()


def cal_fdr (infile):
  '''
  calculate fdr for each of the two populations: pval and pval_b
  '''
  df = pd.read_csv(infile, sep='\t', low_memory=False)

  pvals = df['pval'].tolist()
  fdrs = mu.fdrcorrection(pvals)
  FDR, ACCEPT = fdrs[1], fdrs[0]
  df['pval_fdr'], df['pval_accept'] = FDR, ACCEPT
  
  pvals = df['pval_b'].tolist()
  fdrs = mu.fdrcorrection(pvals)
  FDR, ACCEPT = fdrs[1], fdrs[0]
  df['pvalb_fdr'], df['pvalb_accept'] = FDR, ACCEPT
    
  df.to_csv(infile, sep = '\t', index = False)

def UnitTest ():
  infile = '../UnitTest_summary/GSM1087987_C0FCSb.library_summary.tsv'
  cal_fdr (infile)
  
#UnitTest ()
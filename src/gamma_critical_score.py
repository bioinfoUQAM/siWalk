'''
Date: 2024-04-18
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gamma.html
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html#scipy.stats.rv_continuous.fit

# The "percent point function (ppf)" outputs the value of a random variable such that its probability is less than or equal to an input probability value.
# ppf(q, alpha, mu, beta)
# @q: critical score, z-score or quantile
'''
from scipy.stats import gamma
import pandas as pd

#= global variables
to_test_upper = 'Howell, Howellb, Guo, Guo_b'.split(', ')
to_test_lower = 'pval, pval_b'.split(', ')

def determine_cdf (infile, threshold=0.05):
  ''' Called by summerize_contigs.py '''
  outfile = infile[:-3] + 'cdf.tsv'
  df = pd.read_csv(infile, sep='\t', low_memory=False).sort_values('CONTIG')
  cols = [] # cols = ['pval_fdr']
  for i in to_test_upper:
    df = evaluating_cdf_upper(i, df)
    cols.append(i + '_ccdf')
  for i in to_test_lower:
    df = evaluating_cdf(i, df)
    cols.append(i + '_cdf')
  df['vote'] = df[cols].apply(lambda row: (row < threshold).sum(), axis=1)
  df.to_csv(outfile, sep='\t', index = False)
  return outfile

def calculate_cdf(x, alpha, loc, scale):
    return gamma.cdf(x, alpha, loc=loc, scale=scale) if pd.notnull(x) else None

def evaluating_cdf(feature, df):
  fit_alpha, fit_loc, fit_beta = gamma.fit(df[feature].dropna().astype(float), method='MLE') ## mitigate errors
  df[feature + '_cdf'] = df[feature].apply(lambda x: calculate_cdf(x, fit_alpha, fit_loc, fit_beta))
  return df

def evaluating_cdf_upper(feature, df):
  ''' The upper CDF is alaso referred to as the Survival Function or the Complementary CDF (CCDF). '''
  fit_alpha, fit_loc, fit_beta = gamma.fit(df[feature].dropna().astype(float), method='MLE') ## mitigate errors
  df[feature + '_ccdf'] = df[feature].apply(lambda x: 1 - calculate_cdf(x, fit_alpha, fit_loc, fit_beta))
  return df

if __name__ == "__main__":
  pass
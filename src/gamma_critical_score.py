'''
Gamma-distribution significance thresholds for phasing scores.

Fits a Gamma distribution to each phasing score column using MLE, then
computes per-row CDF (for p-values, where low = significant) or CCDF
(for Howell/Guo scores, where high = significant). Adds a vote column
counting how many scores pass the threshold.

Author: Chao-Jung Wu
Date:   2024-Apr-18
'''
from scipy.stats import gamma
import pandas as pd

to_test_upper = 'Howell, Howellb, Guo, Guo_b'.split(', ')
to_test_lower = 'pval, pval_b'.split(', ')

def determine_cdf(infile, threshold=0.05):
  """Fit Gamma distributions and append CDF/CCDF significance columns. Called by summarize_contigs."""
  outfile = infile[:-3] + 'cdf.tsv'
  df = pd.read_csv(infile, sep='\t', low_memory=False).sort_values('CONTIG')
  cols = []
  for i in to_test_upper:
    df = evaluating_cdf_upper(i, df)
    cols.append(i + '_ccdf')
  for i in to_test_lower:
    df = evaluating_cdf(i, df)
    cols.append(i + '_cdf')
  df['vote'] = df[cols].apply(lambda row: (row < threshold).sum(), axis=1)
  df.to_csv(outfile, sep='\t', index=False)
  return outfile

def calculate_cdf(x, alpha, loc, scale):
  """Return Gamma CDF value, or None if x is NaN."""
  return gamma.cdf(x, alpha, loc=loc, scale=scale) if pd.notnull(x) else None

def evaluating_cdf(feature, df):
  """Fit Gamma to *feature* and append its CDF column (significant when low)."""
  fit_alpha, fit_loc, fit_beta = gamma.fit(df[feature].dropna().astype(float), method='MLE')
  df[feature + '_cdf'] = df[feature].apply(lambda x: calculate_cdf(x, fit_alpha, fit_loc, fit_beta))
  return df

def evaluating_cdf_upper(feature, df):
  """Fit Gamma to *feature* and append its CCDF column (significant when low; aka survival function)."""
  fit_alpha, fit_loc, fit_beta = gamma.fit(df[feature].dropna().astype(float), method='MLE')
  df[feature + '_ccdf'] = df[feature].apply(lambda x: 1 - calculate_cdf(x, fit_alpha, fit_loc, fit_beta))
  return df

if __name__ == "__main__":
  pass

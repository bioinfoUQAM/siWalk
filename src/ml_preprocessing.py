'''
Feature preprocessing for ML input.

Removes non-ML columns (annotations, unnormalized expression features) from
the contig feature table, converts boolean hit columns, defines the class
label, and reorders columns. Outputs a cleaned *.rm_nonfeat.tsv file.

Author: Chao-Jung Wu
Date:   2024-May-16
'''
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import os

def move_A_to_next_B(df, A, B):
  """Move column *A* to immediately after column *B*."""
  cols = list(df.columns)
  cols.pop(cols.index(A))
  cols.insert(cols.index(B) + 1, A)
  return df[cols]

def move_A_to_before_B(df, A, B):
  """Move column *A* to immediately before column *B*."""
  cols = list(df.columns)
  cols.pop(cols.index(A))
  cols.insert(cols.index(B) - 1, A)
  return df[cols]

def non_features():
  """Return columns that are segment annotations, not ML features."""
  cols_to_remove = []
  cols_to_remove += ['chr', 'pos_of_maxf', 'L_bound', 'R_bound', 'star_seq_ifDominantBoth', 'precursor', 'prefold']
  cols_to_remove += ['pval_accept']
  cols_to_remove += ['pvalb_accept', 'pvalb_fdr', 'Literature']
  cols_to_remove += ['start', 'end', 'strand', 'Phas_Ratio', 'Phas_Score', 'Pvalue']
  cols_to_remove += ['Best_miR', 'anyhit', 'twohit']
  cols_to_remove += ['vote']
  cols_to_remove += ['precursor_200_500', 'prefold_200_500', 'precursor_500_200', 'prefold_500_200']
  return cols_to_remove

def expression_features_not_normalized():
  """Return non-normalized expression columns that should be excluded from ML models."""
  cols_to_remove = []
  cols_to_remove += ['p', 'u', 'U', 'maxf']
  cols_to_remove += ['Wfreq_21', 'Cfreq_21', 'cntgfrq_all', 'total_frq_DicerCall']
  cols_to_remove += ['Howell', 'Howellb', 'Guo', 'Guo_b', 'pval, pval_b']
  return cols_to_remove

def run(infile, retainedtag='retained'):
  """Remove non-feature columns from *infile* and write *.rm_nonfeat.tsv*."""
  outfile = infile[:-3] + 'rm_nonfeat.tsv'
  df = pd.read_csv(infile, sep='\t')

  if 'anyhit' in df.columns:
    df['anyhitBool'] = df['anyhit'] != 'na'
  if 'twohit' in df.columns:
    df['twohitBool'] = df['twohit'] != 'na'

  if 'vote' and 'Literature' in df.columns:
    df[retainedtag] = (df['vote'] >= 2) & (df['Literature'] == True)

  keywords = ["Howell >=", "Howellb >=", "Guo >=", "Guo_b >=", "pval_fdr <=", "pval_b <="]
  cols_to_evaluate = []
  for keyword in keywords:
    cols_to_evaluate += [col for col in df.columns if keyword in col]
  cols_to_evaluate += [col for col in df.columns if "pval <=" in col]

  cols_to_remove = cols_to_evaluate + non_features() + expression_features_not_normalized()
  for col in cols_to_remove:
    if col in df.columns:
      df = df.drop([col], axis=1)

  df = move_A_to_next_B(df, 'pval_fdr', 'pval')
  df = move_A_to_before_B(df, 'BestScore', 'anyhitBool')

  df.to_csv(outfile, sep='\t', index=False)
  return outfile

def test_run(infile):
  os.system('cp ' + infile + ' ../output/')
  filename = os.path.basename(infile)
  new_infile = '../output/' + filename
  outfile = run(new_infile)
  return outfile

if __name__ == "__main__":
  infile = '../UnitTest_feature_definition/GSM1087987_C0FCSb.contig_features.addref.threshold.tsv'
  infile = '../UnitTest_ml/____GSM738727_8676.contig_features.tsv'
  outfile = test_run(infile)
  print('\nsee outfile:', outfile)

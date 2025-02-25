'''
2024-May-16
Chao-Jung Wu


backup codes (my pandas notes):
df['twohitBool'] = df['twohit'].apply(lambda x: False if x == 'na' else True)
print(df.columns[10])
'''
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import os

def move_A_to_next_B (df, A, B):
  ''' Move column A to be next to (after) column B '''
  cols = list(df.columns)
  a_index = cols.index(A)
  b_index = cols.index(B)
  # Remove 'A' from its current position
  cols.pop(a_index)
  # Insert 'A' after 'B'
  cols.insert(b_index + 1, A)
  # Reindex DataFrame with the new column order
  df = df[cols]
  return df

def move_A_to_before_B (df, A, B):
  ''' Move column A to be before column B '''
  cols = list(df.columns)
  a_index = cols.index(A)
  b_index = cols.index(B)
  # Remove 'A' from its current position
  cols.pop(a_index)
  # Insert 'A' after 'B'
  cols.insert(b_index - 1, A)
  # Reindex DataFrame with the new column order
  df = df[cols]
  return df

def non_features():
  ''' These columns are annotations to a segment but not features to be learned by machine learning models. '''
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
  ''' Non-normalized expression features can adversely affect machine learning performance across different libraries. Since normalization issues hinder phasing score calculations, these features should be excluded from the models. '''
  cols_to_remove = []
  cols_to_remove += ['p', 'u', 'U', 'maxf']
  cols_to_remove += ['Wfreq_21', 'Cfreq_21', 'cntgfrq_all', 'total_frq_DicerCall']
  cols_to_remove += ['Howell', 'Howellb', 'Guo', 'Guo_b', 'pval, pval_b']
  return cols_to_remove

def run(infile, retainedtag = 'retained'):
    outfile = infile[:-3] + 'rm_nonfeat.tsv'

    df = pd.read_csv(infile, sep = '\t')#; print(df.head()); cols = df.columns.to_list(); print(cols)

    #= Transform some columns
    if 'anyhit' in df.columns:
      df['anyhitBool'] = df['anyhit'] != 'na'
    if 'twohit' in df.columns:
      df['twohitBool'] = df['twohit'] != 'na'
    
    #= Define class based on some rules. This step makes a strong assumption.
    if 'vote' and 'Literature' in df.columns:
      df[retainedtag] = (df['vote'] >= 2) & (df['Literature'] == True)

    keywords = ["Howell >=", "Howellb >=", "Guo >=", "Guo_b >=", "pval_fdr <=", "pval_b <="]
    cols_to_evaluate = []
    for keyword in keywords:
      cols_to_evaluate += [col for col in df.columns if keyword in col]

    #= Remove some columns
    cols_to_evaluate += [col for col in df.columns if "pval <=" in col]
    cols_to_remove = cols_to_evaluate + non_features()
    cols_to_remove += expression_features_not_normalized()
    for col in cols_to_remove:
      if col in df.columns:
        df = df.drop([col], axis=1)

    #= Rearrange cols
    df = move_A_to_next_B (df, 'pval_fdr', 'pval')
    df = move_A_to_before_B (df, 'BestScore', 'anyhitBool')

    #= Save the cleansed dataset to a new file
    df.to_csv(outfile, sep='\t', index = False)
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
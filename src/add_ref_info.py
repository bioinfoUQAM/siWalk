'''
Date: 2024-Mar-26
Update: 2024-Apr-10
This is post-processing, ie adding annotations to the pipeline output.
I want to label contigs that are known to be phasing clusters according to literature.

Make sure the reference file has a header of column names.
'''
import pandas as pd
# import os

def run (ref, infile):
  outfile = infile[:-3] + 'addref.tsv'
  df = pd.read_csv(infile, sep='\t', low_memory=False).sort_values('CONTIG')
  ref = pd.read_csv(ref, sep='\t')
  ref.rename(columns={'CONTIG_literature_seen': 'CONTIG'}, inplace=True)
  ref['Literature'] = True
  if 'chromosome' in ref.columns:
    ref = ref.drop('chromosome', axis=1)
  df_add_ref = pd.merge(df,ref,on='CONTIG',how='left')
  df_add_ref['Literature'].fillna(False, inplace = True)
  df_add_ref.to_csv(outfile, sep='\t', index = False)
  return outfile

if __name__ == "__main__":
  pass
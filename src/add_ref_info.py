'''
Post-processing module: annotates pipeline output with literature-known phasing clusters.

Merges a reference file of known phasing loci with the contig feature table,
adding a boolean 'Literature' column. The reference file must include a header row.

Author: Chao-Jung Wu
Date:   2024-Mar-26
'''
import pandas as pd

def run(ref, infile):
  """Merge *infile* with reference loci and return the annotated output path."""
  outfile = infile[:-3] + 'addref.tsv'
  df = pd.read_csv(infile, sep='\t', low_memory=False).sort_values('CONTIG')
  ref = pd.read_csv(ref, sep='\t')
  ref.rename(columns={'CONTIG_literature_seen': 'CONTIG'}, inplace=True)
  ref['Literature'] = True
  if 'chromosome' in ref.columns:
    ref = ref.drop('chromosome', axis=1)
  df_add_ref = pd.merge(df, ref, on='CONTIG', how='left')
  df_add_ref['Literature'].fillna(False, inplace=True)
  df_add_ref.to_csv(outfile, sep='\t', index=False)
  return outfile

if __name__ == "__main__":
  pass

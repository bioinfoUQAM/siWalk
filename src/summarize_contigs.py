'''
Contig summarization: soft aggregation of phasing windows into segments.

Reads per-contig phasingstat files, applies pre-filter rules (k > 2, maxf > 3),
and aggregates all windows sharing the same effector position (eff_pos) into a
single segment. Retrieves sequences, folds precursors with RNAfold, and appends
FDR and Gamma-CDF significance columns.

Pre-filter: segments with k ≤ 2 or maxf ≤ 3 are discarded.
Precursors are retrieved from the Watson strand.

Author: Chao-Jung Wu
Date:   2024-Apr-05
'''
import os
import pandas as pd
import fdr
import gamma_critical_score as ga
import retrieve_seq as ret


def describe_summary(infile, rep_input, inBasename):
  """Write a describe() summary of segments with k > 2 to a separate TSV."""
  df = pd.read_csv(infile, sep='\t', low_memory=False)
  df = df[df.k > 2]
  describe_file = rep_input + inBasename + '.library_summary.describe.tsv'
  df.describe().to_csv(describe_file, sep='\t')

def keep_consolidated_potential_positions_of_a_contig(infile, outfile, genome_file, DicerCall):
  '''
  Soft-aggregate windows of one contig into segments (current model).

  Groups all qualifying windows by eff_pos. For each group, takes the
  maximum of most features (minimum for L_bound, pval, pval_b) and
  retrieves effector sequence, precursor sequence, and three precursor
  window variants for RNAfold.
  '''
  words = infile.split('/')[-1].split('.')
  if   len(words) == 5: contig = words[2]
  elif len(words) == 4: contig = words[1]
  else:
    print('filename error in summarize_contigs, filename:', infile)
    return

  fho = open(outfile, 'a')
  df = pd.read_csv(infile, sep='\t')
  if len(df) == 0: return

  df = df[df.k > 2]
  df = df[df.maxf > 3]

  cols  = 'k, n, N, p, u, U, maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq, ext_k, L_bound, R_bound, length, Howell, Howellb, Guo, Guo_b, pval, pval_b, dominant_strand, Wfreq, Cfreq, contig_frq, chromosome, mfe'.split(', ')
  cols2 = 'L_bound, pval, pval_b'.split(', ')

  for eff_pos in df['eff_pos'].unique():
    df_tmp = df.loc[df['eff_pos'] == eff_pos]
    values = []
    for col in cols:
      v = max(df_tmp[col])
      if col in cols2: v = min(df_tmp[col])
      if col == 'dominant_strand':
        v = eff_strand = df_tmp[col].value_counts()[:1].index.tolist()[0]
      if col == 'eff_strand':
        v = eff_strand = df_tmp[col].value_counts()[:1].index.tolist()[0]
      if col == 'chromosome':
        v = CHR = df_tmp[col].unique().tolist()[0]
      if col == 'L_bound': start = v
      if col == 'R_bound': end   = v
      values.append(v)
      if col == 'mfe': values.append(min(df_tmp[col]))

    _, eff_seq = ret.retrieve(CHR, eff_pos, eff_pos + DicerCall - 1, genome_file, eff_strand)
    if eff_seq == 'error': continue
    a, b = max(eff_pos - 200, 0), max(eff_pos - 500, 0)
    _, precursor_200_500 = ret.retrieve(CHR, a, eff_pos + 500 - 1, genome_file, eff_strand)
    if precursor_200_500 == 'error': continue
    _, precursor_500_200 = ret.retrieve(CHR, b, eff_pos + 200 - 1, genome_file, eff_strand)
    if precursor_500_200 == 'error': continue
    prefold_200_500, premfe_200_500 = ret.run_RNAfold(precursor_200_500)
    prefold_500_200, premfe_500_200 = ret.run_RNAfold(precursor_500_200)

    star_strand = 'C' if eff_strand == 'W' else 'W'
    star_pos = eff_pos + 2 if eff_strand == 'C' else eff_pos - 2
    _, star_seq_if_dominant_strand_both = ret.retrieve(CHR, star_pos, star_pos + DicerCall - 1, genome_file, star_strand)
    if star_seq_if_dominant_strand_both == 'error': continue

    _, precursor = ret.retrieve(CHR, start, end, genome_file, eff_strand)
    if precursor == 'error': continue
    prefold, premfe = ret.run_RNAfold(precursor)
    dist_5p = precursor.find(eff_seq)
    dist_3p = len(precursor) - int(dist_5p) - len(eff_seq)

    segment = str(CHR) + ':' + str(start) + '-' + str(end)
    mylist  = [contig] + values + [eff_seq, star_seq_if_dominant_strand_both, segment, precursor, prefold, premfe, dist_5p, dist_3p]
    mylist += [precursor_200_500, prefold_200_500, premfe_200_500, precursor_500_200, prefold_500_200, premfe_500_200]
    print('\t'.join([str(x) for x in mylist]), file=fho)

  fho.close()

def run_summary(rep_input, output_tmp, inBasename, genome_file, DicerCall):
  """Aggregate all phasingstat files for one library into a summary TSV with FDR and CDF columns."""
  outfile = output_tmp + inBasename + '.library_summary.tsv'
  with open(outfile, 'w') as fho:
    my_col_names  = 'CONTIG, k, n, N, p, u, U, maxf, pos_of_maxf'
    my_col_names += ', eff_strand, eff_pos, eff_frq'
    my_col_names += ', ext_k, L_bound, R_bound, length'
    my_col_names += ', Howell, Howellb, Guo, Guo_b, pval, pval_b'
    my_col_names += ', dominant_strand, Wfreq_21, Cfreq_21, cntgfrq_all, chr'
    my_col_names += ', max_mfe, min_mfe, eff_seq, star_seq_ifDominantBoth, segment, precursor, prefold, premfe, dist_5p, dist_3p'
    my_col_names += ', precursor_200_500, prefold_200_500, premfe_200_500, precursor_500_200, prefold_500_200, premfe_500_200'
    print('\t'.join([str(x) for x in my_col_names.split(', ')]), file=fho)

  infiles = [f for f in os.listdir(rep_input)
             if f.endswith('.phasingstat.tsv') and os.path.isfile(os.path.join(rep_input, f))]
  for filename in sorted(infiles):
    if inBasename not in filename: continue
    infile = rep_input + filename
    keep_consolidated_potential_positions_of_a_contig(infile, outfile, genome_file, DicerCall)

  fdr.cal_fdr(outfile)
  outfile = ga.determine_cdf(outfile)
  return outfile

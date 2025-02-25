'''
Date: 2024-Apr-05
Extract one meaningful segment from each contig of a library, store in file: x.library_summary.tsv
Then on the summary file, apply describe(), store in file: x.library_summary.describe.tsv

#= Run from my desktop, or
#= Run from Narval:
module load python/3.11.5 # or any python3
module load scipy-stack #= to have pandas

Pre-filter rules:
keep df.k > 2
keep df.maxf > 3
Then the representative segmenet is between min(L_bound) and max(R_bound) of the remaining reporting positions.
With the pre-filter parameters, all the 28 known contigs are present; therefore, if the simple parameters are more stringent, we could lose the known contigs.

Precursors are retrieved always from the Watson strand.
'''
import os
import pandas as pd
import fdr
import gamma_critical_score as ga
import retrieve_seq as ret

# def summary_each_contig (infile, outfile):
  # ''' Hard aggregate (Not Using for the moment)
  # Strickly generate one represetative segment for one contig.
  # This function has not been maintained for a while, it might not work, especially the features need update.
  # '''
  # _, contig, _, _ = infile.split('/')[-1].split('.') # GSM1087987_C0FCSb.1__2785000_2790250.phasingstat.tsv
  # fho = open (outfile, 'a')
  # df = pd.read_csv(infile, sep='\t')
  # if len(df) == 0: return

  # df = df[df.k > 2]
  # df = df[df.maxf > 3]

  # L_bound = min(df.L_bound)
  # R_bound = min(df.R_bound)
  # contig_frq = max(df.contig_frq)
  # freq = max(df.freq)
  # k = max(df.k)
  # n = max(df.n)
  # N = max(df.N)
  # p = max(df.p)
  # u = max(df.u)
  # U = max(df.U)
  # maxf = max(df.maxf)
  # pos_of_maxf = df.loc[df.maxf == maxf, 'pos_of_maxf'].value_counts()[:1].index.tolist()[0]
  # eff_frq = max(df.eff_frq)
  # eff_pos = df.loc[df.eff_frq == eff_frq, 'eff_pos'].value_counts()[:1].index.tolist()[0]
  # eff_strand = df.loc[df.eff_frq == eff_frq, 'eff_strand'].value_counts()[:1].index.tolist()[0]
  # eff_seq = df.loc[df.eff_frq == eff_frq, 'eff_seq'].value_counts()[:1].index.tolist()[0]
  # ext_k = max(df.ext_k)
  # length = R_bound - L_bound + 1
  # Wfreq = max(df.Wfreq)
  # Cfreq = max(df.Cfreq)
  # chromosome = df.chromosome.value_counts()[:1].index.tolist()[0]
  # Howell = max(df.Howell)
  # Howellb = max(df.Howellb)
  # Guo = max(df.Guo)
  # Guo_b = max(df.Guo_b)
  # pval = min(df.pval)
  # pval_b = min(df.pval_b)
  # strand = df['dominant_strand'].value_counts()[:1].index.tolist()[0] # The highest count of the strand is returned.

  # mylist = [contig, freq, k, n, N, p, u, U, maxf, pos_of_maxf,
               # eff_strand, eff_pos, eff_frq, eff_seq, 
               # ext_k, L_bound, R_bound, length,
               # Howell, Howellb, Guo, Guo_b, pval, pval_b,
               # strand, Wfreq, Cfreq, contig_frq, chromosome]
  # print('\t'.join([str(x) for x in mylist]), file=fho)
  # fho.close()
  # return

def describe_summary (infile, rep_input, inBasename):
  df = pd.read_csv(infile, sep='\t', low_memory=False)
  df = df[df.k > 2]
  describe_file = rep_input + inBasename + '.library_summary.describe.tsv'
  df.describe().to_csv(describe_file, sep = '\t')

def keep_consolidated_potential_positions_of_a_contig (infile, outfile, genome_file, DicerCall):
  ''' Soft aggregate (Current using model), Date: 2024-Apr-13
  Collects all segments of one contig if the segment fulfills some criteria.
  Then, segments with the same "active siRNA position (eff_pos)" are aggregated into a single segment, with the boundaries defined by the leftmost and rightmost coordinates. '''
  words = infile.split('/')[-1].split('.')
  if len(words) == 5: contig = words[2]   # S002B4B.6D.1__2785000_2790250.phasingstat.tsv
  elif len(words) == 4: contig = words[1] # GSM1087987_C0FCSb.1__2785000_2790250.phasingstat.tsv
  else:
      print('filename error in summarize_contigs, filename:', infile)
      return

  fho = open (outfile, 'a')
  df = pd.read_csv(infile, sep='\t')
  if len(df) == 0: return

  df = df[df.k > 2]
  df = df[df.maxf > 3]

  cols = 'k, n, N, p, u, U, maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq, ext_k, L_bound, R_bound, length, Howell, Howellb, Guo, Guo_b, pval, pval_b, dominant_strand, Wfreq, Cfreq, contig_frq, chromosome, mfe'.split(', ')

  cols2 = 'L_bound, pval, pval_b'.split(', ')

  eff_pos_values = df['eff_pos'].unique()
  for eff_pos in eff_pos_values:
    df_tmp = df.loc[df['eff_pos'] == eff_pos]
    values = []
    for col in cols:
      v = max(df_tmp[col])
      if col in cols2: v = min(df_tmp[col])
      if col == 'dominant_strand':
        v = dominant_strand = df_tmp[col].value_counts()[:1].index.tolist()[0] # The highest count
      if col == 'eff_strand':
        v = eff_strand = df_tmp[col].value_counts()[:1].index.tolist()[0] # The highest count
      if col == 'chromosome':
        v = CHR = df_tmp[col].unique().tolist()[0] # copy as is
      if col == 'L_bound': start = v
      if col == 'R_bound': end = v
      values.append(v)
      if col == 'mfe': values.append(min(df_tmp[col]))

    _, eff_seq = ret.retrieve (CHR, eff_pos, eff_pos+DicerCall-1, genome_file, eff_strand)
    if eff_seq == 'error': continue
    a, b = max (eff_pos-200, 0), max (eff_pos-500, 0) # ensure left coordinates >= 0
    _, precursor_200_500 = ret.retrieve (CHR, a, eff_pos+500-1, genome_file, eff_strand)
    if precursor_200_500 == 'error': continue
    _, precursor_500_200 = ret.retrieve (CHR, b, eff_pos+200-1, genome_file, eff_strand)
    if precursor_500_200 == 'error': continue
    prefold_200_500, premfe_200_500 = ret.run_RNAfold (precursor_200_500)
    prefold_500_200, premfe_500_200 = ret.run_RNAfold (precursor_500_200)
    
    star_strand = 'C' if eff_strand == 'W' else 'W'
    star_pos = eff_pos + 2 if eff_strand == 'C' else eff_pos - 2
    _, star_seq_if_dominant_strand_both = ret.retrieve (CHR, star_pos, star_pos+DicerCall-1, genome_file, star_strand)
    if star_seq_if_dominant_strand_both == 'error': continue

    _, precursor = ret.retrieve (CHR, start, end, genome_file, eff_strand)
    if precursor == 'error': continue
    prefold, premfe = ret.run_RNAfold (precursor)
    dist_5p = precursor.find(eff_seq) # location of eff_seq within precursor, relative to 5p
    dist_3p = len(precursor) - int(dist_5p) - len(eff_seq)

    segment = str(CHR) + ':' + str(start) + '-' + str(end)

    mylist = [contig] + values + [eff_seq, star_seq_if_dominant_strand_both, segment, precursor, prefold, premfe, dist_5p, dist_3p]
    mylist += [precursor_200_500, prefold_200_500, premfe_200_500, precursor_500_200, prefold_500_200, premfe_500_200]
    print('\t'.join([str(x) for x in mylist]), file=fho)

  fho.close()
  return

def run_summary (rep_input, output_tmp, inBasename, genome_file, DicerCall):
  outfile = output_tmp + inBasename + '.library_summary.tsv'
  fho = open (outfile, 'w')
  my_col_names = 'CONTIG, k, n, N, p, u, U, maxf, pos_of_maxf'
  my_col_names += ', eff_strand, eff_pos, eff_frq'
  my_col_names += ', ext_k, L_bound, R_bound, length'
  my_col_names += ', Howell, Howellb, Guo, Guo_b, pval, pval_b'
  my_col_names += ', dominant_strand, Wfreq_21, Cfreq_21, cntgfrq_all, chr'
  my_col_names += ', max_mfe, min_mfe, eff_seq, star_seq_ifDominantBoth, segment, precursor, prefold, premfe, dist_5p, dist_3p'
  my_col_names += ', precursor_200_500, prefold_200_500, premfe_200_500, precursor_500_200, prefold_500_200, premfe_500_200'

  print('\t'.join([str(x) for x in my_col_names.split(', ')]), file=fho)
  fho.close()

  infiles = [f for f in os.listdir(rep_input) if f.endswith('.phasingstat.tsv') and os.path.isfile(os.path.join(rep_input, f))]
  for filename in sorted(infiles):
    if inBasename not in filename: continue
    infile = rep_input + filename
    keep_consolidated_potential_positions_of_a_contig (infile, outfile, genome_file, DicerCall)

  fdr.cal_fdr (outfile)
  outfile = ga.determine_cdf(outfile)
  return outfile
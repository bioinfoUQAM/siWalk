'''
Effector siRNA localization prediction within a precursor sequence (Module B).

Given a precursor RNA sequence, exhaustively generates candidate siRNA windows
(DicerCall ± 2 nt), extracts structural features via RNAfold, runs miRanda
target search, applies an ML classifier for position scoring, and selects the
top-ranked start/end position by a correlation-weighted indication score.

Entry points:
- run_one_precursor(): programmatic API returning (start, end, score, outfiles)
- user_interface(): CLI wrapper that renames output files and draws a plot

Author: Chao-Jung Wu
'''
from datetime import datetime
import sys, os
import numpy as np
import random
import pandas as pd
from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder, QuantileTransformer

import retrieve_seq as ret
import create_more_features as cmf
import mirCheck_eval_hairpin as hp
import miRanda_search_target as mst
import barplot_indication as bpi
import mlclassify_localization as mlloc


def set_seed(seed=0):
  """Set Python and NumPy random seeds for reproducibility."""
  random.seed(seed)
  np.random.seed(seed)


seed = 0
set_seed(seed=seed)
tag = 'consistent'

k = 100
file_used_as_training = '../dbs/background.tsv'
datafile = '../model/Arabidopsis_strcture_feature_importance_n_correlation.tsv'

cols_to_drop = ['CONTIG', 'eff_seq', 'retained', tag, 'segment']
cols2drop_for_locapredi = [
'CONTIG', 'k', 'n', 'N', 'eff_strand', 'eff_pos', 'eff_frq', 'ext_k', 'length',
'pval', 'pval_fdr', 'pval_b', 'dominant_strand', 'max_mfe', 'min_mfe', 'eff_seq',
'segment', 'premfe_200_500', 'premfe_500_200', 'Howell_ccdf', 'Howellb_ccdf',
'Guo_ccdf', 'Guo_b_ccdf', 'pval_cdf', 'pval_b_cdf', 'A%', 'C%', 'G%', 'T%', 'GC%',
'5p1A', '5p1C', '5p1G', '5p1T', '5p2A', '5p2C', '5p2G', '5p2T', '5p3A', '5p3C',
'5p3G', '5p3T', '5p4A', '5p4C', '5p4G', '5p4T', '5p5A', '5p5C', '5p5G', '5p5T',
'3p1A', '3p1C', '3p1G', '3p1T', '3p2A', '3p2C', '3p2G', '3p2T', '3p3A', '3p3C',
'3p3G', '3p3T', '3p4A', '3p4C', '3p4G', '3p4T', '3p5A', '3p5C', '3p5G', '3p5T',
'md1A', 'md1C', 'md1G', 'md1T', 'md2A', 'md2C', 'md2G', 'md2T', 'md3A', 'md3C',
'md3G', 'md3T', 'md4A', 'md4C', 'md4G', 'md4T', 'md5A', 'md5C', 'md5G', 'md5T',
'mircheck_conclu25', 'fback_start25', 'fback_stop25', 'mircheck_conclu52',
'fback_start52', 'fback_stop52', 'retained', tag]
cols_to_drop += cols2drop_for_locapredi


def score(p, l, L):
  """Look up the weighted indication score for start position *p* and length *l* in *L*."""
  for i in L:
    if int(i[0]) == p and int(i[1]) == l:
      return float(i[2])
  return 0


def S_score(p, L):
  """Compute the start-position score S(p) summed over DicerCall-range lengths."""
  Sp = 0
  for l in range(19, 24):
    Sp += score(p, l, L)
  return Sp


def E_score(p, L):
  """Compute the end-position score E(p) summed over DicerCall-range lengths."""
  Ep = 0
  for l in [19, 20, 21, 22, 23]:
    Ep += score(p - l + 1, l, L)
  return Ep


def argmax_local(p, L):
  """Return (max_score, best_length) for position *p* across DicerCall-range lengths."""
  maximun, best_longeur = -10000, 0
  for l in range(19, 24):
    value = S_score(p, L) + E_score(p + l - 1, L)
    if value > maximun:
      maximun = value
      best_longeur = l
  return maximun, best_longeur


def argmax_global(pspace, L):
  """Score all positions in *pspace* and return rows sorted by localMax descending."""
  data = []
  for p in pspace:
    Sp = S_score(p, L)
    Ep = E_score(p, L)
    localMax, local_best_longeur = argmax_local(p, L)
    End_position = p + local_best_longeur - 1
    data.append([p, Sp, Ep, localMax, local_best_longeur, End_position])
  data.sort(reverse=True, key=lambda x: x[3])
  return data


def compute_indications_for_effector_start_end(infile, position_list):
  """Select the best effector start/end from an indication score file.

  Returns (start, end, globalMax, full_indication_tsv, top6_recommendation_tsv).
  ML-predicted positions are prioritised in the top-6 recommendation.
  """
  outfile = infile[:-3] + 'indication.tsv'
  with open(infile, 'r') as fh:
    L = [x.rstrip('\n').split('\t') for x in fh.readlines()][1:]
  pspace = sorted(list(set([int(i[0]) for i in L])))
  data = argmax_global(pspace, L)

  for item in data:
    item.append(item[0] in position_list)

  fho = open(outfile, 'w')
  print('\t'.join('Position p, Start S(p), End E(p), Sum of Indications, Best length, End position, ML predicted'.split(', ')), file=fho)
  for i in data:
    print('\t'.join([str(round(x, 6)) for x in i]), file=fho)
  fho.close()

  if len(position_list) > 0:
    filtered = [item for item in data if item[0] in position_list]
    remaining = [item for item in data if item[0] not in position_list]
    data = filtered + remaining
    data = data[:6]
  outfile2 = infile[:-3] + 'top6_recommendation.tsv'
  fho = open(outfile2, 'w')
  print('\t'.join('Position p, Start S(p), End E(p), Sum of Indications, Best length, End position, ML predicted'.split(', ')), file=fho)
  for i in data:
    print('\t'.join([str(round(x, 6)) for x in i]), file=fho)
  fho.close()

  start, _, _, globalMax, _, end, _ = data[:1][0]
  return start, end, round(globalMax, 6), outfile, outfile2


def pretreat_location_features(infile_of_Instance):
  """Prepend placeholder label columns and concatenate with the training background TSV."""
  outfile = infile_of_Instance[:-3] + 'pretreated.tsv'
  df1 = pd.read_csv(infile_of_Instance, sep='\t', index_col=False)
  NB_instances = len(df1)
  df1['CONTIG'] = 'nana'
  df1['retained'] = True
  df1['consistent'] = True
  df2 = pd.read_csv(file_used_as_training, sep='\t', index_col=False)
  df2['_longeur_'] = 0
  result_df = concatenate(df1, df2)
  result_df.to_csv(outfile, sep='\t', index=False)
  return outfile, NB_instances


def concatenate(df1, df2):
  """Concatenate df1 and df2, aligning to df2's column order and filling missing values with 0.

  Parameters:
  df1, potentially with fewer columns.
  df2, with the full set of columns.

  Returns:
  The concatenated DataFrame with the same columns and order as df2, with missing values filled.

  Example:
  >>> df1 = pd.DataFrame({
  ...     'A': [1],
  ...     'B': [4]
  ... })
  >>> df2 = pd.DataFrame({
  ...     'A': [2, 3],
  ...     'B': [5, 6],
  ...     'C': [7, 8],
  ...     'D': [9, 10]
  ... })
  >>> cc(df1, df2)
     A  B    C     D
  0  1  4  0.0  0.0
  1  2  5  7.0  9.0
  2  3  6  8.0 10.0
  """
  result_df = pd.concat([df1, df2], ignore_index=True)
  result_df = result_df.reindex(columns=df2.columns)
  result_df = result_df.fillna(0)
  return result_df


def encode_and_compute_weight(infile, NB_instances, datafile, k=k, cols_to_drop=cols_to_drop):
  """Encode structural features and compute correlation-weighted indication scores."""
  output_file = infile[:-3] + 'weight.tsv'
  feature_names_topk, feature_correlations, feature_importances = get_data(datafile, k)

  df = pd.read_csv(infile, sep='\t', index_col=False)
  Xs = df.drop(columns=cols_to_drop, axis=1, errors='ignore')

  ord_col = []
  cat_col = [i for i, var in enumerate(Xs.columns) if Xs[var].dtype == 'O' or Xs[var].dtype == 'bool']
  num_col = [i for i, var in enumerate(Xs.columns) if Xs[var].dtype != 'O' and Xs[var].dtype != 'bool']
  Xs_encoded = preprocess_direct(Xs, num_col, cat_col, ord_col)
  Xs_encoded = Xs_encoded.head(NB_instances)
  df = df.head(NB_instances)

  weighted_sum_localpredi = calculate_weighted_sum_based_on_correlation(Xs_encoded, feature_importances, feature_correlations)
  df['weighted_sum_localpredi'] = weighted_sum_localpredi

  df = df[['dist_5p', '_longeur_', 'weighted_sum_localpredi']]
  df.to_csv(output_file, sep='\t', index=False)
  return output_file


def calculate_weighted_sum_based_on_correlation(df, feature_importances, feature_correlations):
  """Calculate the correlation × importance weighted sum of features for each instance.

  Adding 0.1 to feature values ensures the product term never becomes zero for
  features with values at or near zero, boosting their contribution to the score.

  Returns an np.array of weighted sums, one per instance.
  """
  weighted_sum = np.zeros(len(df))
  for feature, importance in feature_importances.items():
    if feature in df.columns and feature in feature_correlations:
      correlation = feature_correlations[feature]
      weighted_sum += correlation * importance * df[feature]
  return weighted_sum


def get_top_k_features(chi2_scores_dict, k=100):
  """Return the top-k features sorted by chi2 score (descending)."""
  sorted_features = sorted(chi2_scores_dict.items(), key=lambda item: item[1], reverse=True)
  top_k_features = dict(sorted_features[:k])
  return top_k_features


def preprocess_direct(X, num_col, cat_col, ord_col):
  """Encode and scale features: quantile-transform numerical, one-hot categorical."""
  if all(isinstance(col, int) for col in num_col):
    num_col = X.columns[num_col]
  if all(isinstance(col, int) for col in cat_col):
    cat_col = X.columns[cat_col]
  if all(isinstance(col, int) for col in ord_col):
    ord_col = X.columns[ord_col]

  scaler = QuantileTransformer(n_quantiles=10, random_state=seed)
  cat_encode = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
  ord_encode = OrdinalEncoder()

  X[num_col] = scaler.fit_transform(X[num_col])
  encoded_cat_array = cat_encode.fit_transform(X[cat_col])
  encoded_cat_df = pd.DataFrame(encoded_cat_array, columns=cat_encode.get_feature_names_out(cat_col), index=X.index)

  X = X.drop(columns=cat_col)
  X = pd.concat([X, encoded_cat_df], axis=1)
  X[ord_col] = ord_encode.fit_transform(X[ord_col])
  return X


def get_data(datafile, k=1000):
  """Load feature importances and correlations from a feature evaluation TSV."""
  df = pd.read_csv(datafile, sep='\t')
  feature_importances = dict(zip(df['feature'], df['importance(chi2)']))
  feature_correlations = dict(zip(df['feature'], df['correlation(pointbiserialr)']))
  feature_names_topk = get_top_k_features(feature_importances, k).keys()
  return feature_names_topk, feature_correlations, feature_importances


def miRNA_target_search(precursor, output_tmp='../tmp/', species='ath', mirbase_file='../dbs/mature.fa'):
  """Run miRanda target search on *precursor* and return (Best_miR, BestScore, anyhit, twohit)."""
  if not os.path.exists(output_tmp): os.makedirs(output_tmp)
  timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
  precursor_name = 'precursorName'
  xfile = output_tmp + precursor_name + '_' + timestamp + '.tsv'
  e = ['mySegmentName', precursor]
  with open(xfile, 'w') as fh:
    print('segment\tprecursor', file=fh)
    print('\t'.join(e), file=fh)
  mRnd_obj = mst.miRanda_class(xfile, mirbase_file, species, output_tmp, output_tmp, precursor_name)
  _, _, Best_miR, BestScore, anyhit_data, twohit_data = mRnd_obj.search_trigger(e)
  return Best_miR, BestScore, anyhit_data, twohit_data


def create_structure_features(precursor, seq2test_effseq, common):
  """Generate structural features for a candidate siRNA window within *precursor*.

  Parameters:
  -----------
  precursor : str
      The full RNA sequence from which the effective sequence is derived.
  seq2test_effseq : str
      The effective sequence (a substring of the precursor) for which
      structural features are to be calculated.

  Returns:
  --------
  df : A DataFrame containing the calculated structural features.
  """
  prefold, premfe, Best_miR, BestScore, anyhit_data, twohit_data = common

  eff_seq = seq2test_effseq
  df = pd.DataFrame({'eff_seq': [eff_seq]})
  df = cmf.mers123(df)

  data = {}
  data['premfe'] = premfe

  start = dist_5p = precursor.find(eff_seq)
  stop = dist_5p + len(eff_seq)
  dist_3p = len(precursor) - int(dist_5p) - len(eff_seq)
  data['dist_5p'], data['dist_3p'] = dist_5p, dist_3p

  a = max(start - 4, 0)
  b = min(stop + 4, len(precursor))
  substring_fold = prefold[a: b]
  ref_seq = precursor[a: b]

  data['paired_percentage'] = cmf.get_paired_percentage(substring_fold)
  data['length_longest_bulge'] = cmf.length_largest_bulge(substring_fold)
  data['length_longest_loop'] = cmf.length_longest_loop(substring_fold)
  data['longest_paired_length'] = cmf.length_largest_bracket_sequence(substring_fold)
  data['mircheck_conclu'], data['fback_start'], data['fback_stop'] = hp.call_mirCheck(prefold, start, stop)

  for window in [3, 5, 7]:
    data['paired_roll' + str(window)] = cmf.get_paired_rolling_average(substring_fold, window)

  for motif in ['AAA', 'TTT', 'CCC', 'GGG']:
    data['NBtriplet' + motif[0]] = cmf.number_of_motif(ref_seq, motif)

  ds = cmf.main_outtsv(ref_seq, substring_fold)
  for d in ds:
    for k, v in d.items():
      if start < 0:
        data[k] = -1
      else:
        data[k] = v

  data['BestScore'] = BestScore
  data['anyhitBool'] = anyhit_data != 'na'
  data['twohitBool'] = twohit_data != 'na'

  new_df = pd.DataFrame([data])
  df = df.drop(['eff_seq'], axis=1)
  df = pd.concat([df, new_df], axis=1)
  return df


def get_siRNA_structure(name, precursor, DicerCall=21, tmpdir='../tmp/', species='ath', mirbase_file='../dbs/mature.fa'):
  """Generate structural features for all candidate siRNA windows in *precursor*.

  Iterates over all start positions and DicerCall ± 2 lengths, builds a feature
  DataFrame, and writes it to a timestamped TSV. Returns the TSV path.
  """
  if not os.path.exists(tmpdir): os.makedirs(tmpdir)
  timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

  prefold, premfe = ret.run_RNAfold(precursor)
  Best_miR, BestScore, anyhit_data, twohit_data = miRNA_target_search(precursor, tmpdir, species, mirbase_file)
  common = [prefold, premfe, Best_miR, BestScore, anyhit_data, twohit_data]

  df = pd.DataFrame()
  for i in range(len(precursor)-DicerCall+1):
    for j in range(DicerCall-2, DicerCall+3):  # DicerCall ± 2 drift range: [19, 20, 21, 22, 23]
      sseq = precursor[i:i+j]
      df_tmp = create_structure_features(precursor, sseq, common)
      df_tmp['_longeur_'] = j
      df = pd.concat([df, df_tmp])
  siRNA_structure_file = tmpdir + timestamp + name + '_siRNA_structures.tsv'
  df.to_csv(siRNA_structure_file, sep='\t', index=False)
  return siRNA_structure_file


def run_one_precursor(name, precursor, DicerCall=21,
                      tmpdir='../tmp/', species='ath', mirbase_file='../dbs/mature.fa'):
  """Run the full localization prediction pipeline for one precursor sequence.

  Returns (start, end, score, indication_tsv, top6_recommendation_tsv).
  """
  siRNA_structure_file = get_siRNA_structure(name, precursor, DicerCall, tmpdir, species, mirbase_file)
  position_list = mlloc.classify_a_file(siRNA_structure_file)
  te, NB_instances = pretreat_location_features(siRNA_structure_file)
  siRNA_indication_file = encode_and_compute_weight(te, NB_instances, datafile)
  start, end, score, outfile, outfile2 = compute_indications_for_effector_start_end(siRNA_indication_file, position_list)
  print(f"{name}, predicted siRNA start: {start}, end: {end}, score: {score}")
  return start, end, score, outfile, outfile2


def dna_to_rna(dna_sequence):
  """Convert a DNA sequence to RNA by replacing T with U."""
  rna_sequence = dna_sequence.replace('T', 'U').replace('t', 'u')
  return rna_sequence


def rna_to_dna(rna_sequence):
  """Convert an RNA sequence to DNA by replacing U with T."""
  dna_sequence = rna_sequence.replace('U', 'T').replace('u', 't')
  return dna_sequence


def user_interface(name, pri, DicerCall, outdir):
  """Run localization prediction and save renamed output files with a candidate plot."""
  start, end, score, outfile, outfile2 = run_one_precursor(name, pri, DicerCall)
  outfile_newname = outdir + name + '.effector_localization_indication.tsv'
  outfile2_newname = outdir + name + '.effector_localization_top6_recommendation.tsv'
  os.rename(outfile, outfile_newname)
  os.rename(outfile2, outfile2_newname)
  bpi.draw_6candidates_interface(outfile_newname, outfile2_newname, pri)


if __name__ == "__main__":
  timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
  outdir = '../output/' + timestamp + '/'; os.makedirs(outdir)
  precursorName = 'yourPrecursor'
  args = sys.argv[1:]
  if len(args) > 2 or len(args) == 0:
    print("The script take two arguments: the sequence of siRNA generating locus; and DicerCall")
    print("usage: python siWalk_predict_siRNA_location.py $priseq $DicerCall")
    print("example usage: TAS3 (PHAS21-21) segment 3__5862036_5862355 with siRNA=TTCTTGACCTTGTAAGACCCC located between 50 and 70 ")
    print("priseq=TCTAGATGATGCATTTCATTATTCTCTTTTTCTTGACCTTGTAAGGCCTTTTCTTGACCTTGTAAGACCCCATCTCTTTCTAAACGTTTTATTATTTTCTCGTTTTACAGATTCTATTCTA")
    print("example usage: TAS3 (PHAS21-21) segment 3_5862187_5862334 with siRNA=TTCTTGACCTTGTAAGACCCC located between 19 and 39 ")
    print("priseq=CTTGACCTTGTAAGGCCTTTTCTTGACCTTGTAAGACCCCATCTCTTTCTAAACGTTTTATTATTTTCTCGTTTTACAGATTCTATTCTATCTCTTCTCAATATAGAATAGATATCTATCT")
    print("DicerCall=21")
    print("python siWalk_predict_siRNA_location.py $priseq $DicerCall")
    sys.exit(0)

  pri, DicerCall = args[0], int(args[1])
  pri = rna_to_dna(pri)  # convert RNA to DNA if needed
  user_interface(precursorName, pri, DicerCall, outdir)
  print('==== See results in', outdir)
  pass

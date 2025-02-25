'''
Chao-Jung Wu
2024-08-13

If need to run on compute canada:

module load StdEnv/2020 viennarna/2.5.1
module load gcc/9.3.0 blast+/2.14.0
module load samtools/1.17 bowtie/1.3.0
module load scipy-stack
#virtualenv --no-download env
source env/bin/activate
#pip install --no-index --upgrade pip
#pip install --no-index statsmodels
#pip install --no-index sklearn
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
  ''' random comes from many places. Here my pipeline could be affected by python and array randomness. Set seed to make the experiments reproducible. '''
  random.seed(seed) #python's randomness
  np.random.seed(seed) # random in arrays

seed = 0
set_seed(seed=seed)
tag='consistent'

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



def score (p, l, L):
    ''' @L: the list of scores (ie. weighted_sum_localpredi = sigma {correlation * importance * df.feature}) '''
    for i in L:
      if int(i[0]) == p and int(i[1]) == l:
          return float(i[2])
    return 0

def S_score(p, L):
    Sp = 0
    for l in range(19, 24):
        Sp += score(p, l, L)
    return Sp

def E_score(p, L):
    Ep = 0
    for l in [19, 20, 21, 22, 23]:
        Ep += score(p - l + 1, l, L)
    return Ep

def argmax_local (p, L):
  maximun, best_longeur = -10000, 0
  for l in range(19, 24):
      value = S_score(p, L) + E_score(p + l -1, L)
      if value > maximun:
          maximun = value
          best_longeur = l
  return maximun, best_longeur

def argmax_global(pspace, L):
    data = []
    for p in pspace:
        Sp = S_score(p, L)
        Ep = E_score(p, L)
        localMax, local_best_longeur = argmax_local (p, L)
        End_position = p+local_best_longeur-1
        data.append([p, Sp, Ep, localMax, local_best_longeur, End_position])
    data.sort(reverse=True, key=lambda x: x[3])  # Sort by localMax from high to low
    return data

def compute_indications_for_effector_start_end(infile, position_list):
    ''' @return: the most likely start and end positions of effector on the given precursor 
    data items  [p,  Sp,     Ep,      localMax, local_best_longeur, End_position]
    data values [38, 356.29, -100.95, 749.86,   22,                 59]
    '''
    outfile = infile[:-3] + 'indication.tsv'
    with open (infile, 'r') as fh:
        L = [x.rstrip('\n').split('\t') for x in fh.readlines()][1:]
    pspace = sorted(list(set([int(i[0]) for i in L])))
    data = argmax_global(pspace, L)  

    # Append "ML predicted" to each item in data
    for item in data:
      item.append(item[0] in position_list)

    fho = open (outfile, 'w')
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
    fho = open (outfile2, 'w')
    print('\t'.join('Position p, Start S(p), End E(p), Sum of Indications, Best length, End position, ML predicted'.split(', ')), file=fho)
    for i in data:
        print('\t'.join([str(round(x, 6)) for x in i]), file=fho)
    fho.close()

    start, _, _,globalMax,_,end,_ = data[:1][0]
    return start, end, round(globalMax, 6), outfile, outfile2

def pretreat_location_features(infile_of_Instance):
    outfile = infile_of_Instance[:-3] + 'pretreated.tsv'
    df1 = pd.read_csv(infile_of_Instance, sep = '\t', index_col=False)
    NB_instances = len(df1)
    df1['CONTIG'] = 'nana'
    df1['retained'] = True
    df1['consistent'] = True
    df2 = pd.read_csv(file_used_as_training, sep = '\t', index_col=False)
    df2['_longeur_'] = 0
    result_df = concatenate(df1, df2)
    result_df.to_csv(outfile, sep='\t', index=False)
    return outfile, NB_instances

def concatenate(df1, df2):
    """ Concatenate two DataFrames, df1 and df2, while maintaining the same columns 
    and order as in df2. If df1 has fewer features (a subset of df2), the function 
    ensures that the resulting DataFrame includes all columns from df2, with missing 
    values filled as needed.

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
    # Reindex to ensure it has the same columns as df2
    result_df = result_df.reindex(columns=df2.columns)
    # Fill any missing values if needed (e.g., NaN values in columns from df2 that were not in df1)
    result_df = result_df.fillna(0)  # or any other value
    return result_df

def encode_and_compute_weight(infile, NB_instances, datafile, k=k, cols_to_drop=cols_to_drop):
  output_file = infile[:-3] + 'weight.tsv'
  feature_names_topk, feature_correlations, feature_importances = get_data(datafile, k)  

  df = pd.read_csv(infile, sep = '\t', index_col=False)
  # Isolate structural features for siRNA localization prediction, excluding expressional ones; Ignore any missing columns
  Xs = df.drop(columns=cols_to_drop,axis=1, errors='ignore')

  ord_col = []
  cat_col = [i for i, var in enumerate(Xs.columns) if Xs[var].dtype == 'O' or Xs[var].dtype == 'bool']
  num_col = [i for i, var in enumerate(Xs.columns) if Xs[var].dtype != 'O' and Xs[var].dtype != 'bool']
  Xs_encoded = preprocess_direct(Xs, num_col, cat_col, ord_col)
  Xs_encoded = Xs_encoded.head(NB_instances)
  df = df.head(NB_instances)

  weighted_sum_localpredi = calculate_weighted_sum_based_on_correlation(Xs_encoded, feature_importances, feature_correlations)
  df['weighted_sum_localpredi'] = weighted_sum_localpredi

  # Clean the DataFrame and retain the relevant information 
  df = df[['dist_5p', '_longeur_', 'weighted_sum_localpredi']] 
  df.to_csv(output_file, sep='\t', index=False)
  return output_file

# def run_indication(siRNA_structure_file):
  # te, NB_instances = pretreat_location_features(siRNA_structure_file)
  # siRNA_indication_file = encode_and_compute_weight(te, NB_instances, datafile)
  # start, end, score, outfile = compute_indications_for_effector_start_end(siRNA_indication_file)
  # return start, end, score, outfile

#######################################################
#######################################################
def calculate_weighted_sum_based_on_correlation(df, feature_importances, feature_correlations):
    """
    Calculate the weighted sum of features based on their importances and correlations with the class.
    Args:
        df (pd.DataFrame): The dataframe containing the features.
        feature_importances (dict): A dictionary where keys are feature names and values are their importances.
        feature_correlations (dict): A dictionary where keys are feature names and values are correlations.
    Returns:
        weighted_sum (np.array): Array of weighted sums for each instance. 

    Effect of "df[feature]+0.1": If the values in df[feature] are small or close to zero, adding 0.1 ensures that the product term (df[feature] + 0.1) never becomes zero (df[feature] values are greater or equal to zero). This boosts all the values slightly and can especially affect features with smaller values, increasing their contribution to the weighted sum.
    """
    weighted_sum = np.zeros(len(df))  
    for feature, importance in feature_importances.items():
        if feature in df.columns and feature in feature_correlations:
            correlation = feature_correlations[feature]
            weighted_sum += correlation * importance * df[feature]
    return weighted_sum

def get_top_k_features(chi2_scores_dict, k=100):
    # Sort the dictionary by chi2 scores in descending order and extract the top k features
    sorted_features = sorted(chi2_scores_dict.items(), key=lambda item: item[1], reverse=True)
    top_k_features = dict(sorted_features[:k])
    return top_k_features

def preprocess_direct(X, num_col, cat_col, ord_col):
    # Check if column lists are indices or names and convert indices to names
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
   df = pd.read_csv(datafile, sep='\t')
   feature_importances = dict(zip(df['feature'], df['importance(chi2)'])) #importance(mi)
   feature_correlations = dict(zip(df['feature'], df['correlation(pointbiserialr)']))
   feature_names_topk = get_top_k_features(feature_importances, k).keys()
   return feature_names_topk, feature_correlations, feature_importances

def miRNA_target_search(precursor, output_tmp = '../tmp/', species='ath', mirbase_file='../dbs/mature.fa'):
    if not os.path.exists(output_tmp): os.makedirs(output_tmp)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    precursor_name = 'precursorName'
    xfile = output_tmp + precursor_name + '_'+ timestamp +'.tsv'
    e = ['mySegmentName', precursor]
    with open(xfile, 'w') as fh:
        print('segment\tprecursor', file=fh)
        print('\t'.join(e), file=fh)
    mRnd_obj = mst.miRanda_class(xfile, mirbase_file, species, output_tmp, output_tmp, precursor_name)
    _, _, Best_miR, BestScore, anyhit_data, twohit_data = mRnd_obj.search_trigger(e)
    return Best_miR, BestScore, anyhit_data, twohit_data

def create_structure_features(precursor, seq2test_effseq, common): #, output_tmp = '../tmp/', species='ath', mirbase_file='../dbs/mature.fa'):
    '''
    Given a precursor sequence and a substring of that precursor, this function 
    generates various structural features related to RNA folding and motif occurrences. 
    It returns these features in a DataFrame.

    Parameters:
    -----------
    precursor : str
        The full RNA sequence from which the effective sequence is derived.
    seq2test_effseq : str
        The effective sequence (a substring of the precursor) for which the 
        structural features are to be calculated.

    Returns:
    --------
    df : A DataFrame containing the calculated structural features of the 
        effective sequence within the precursor.

    Example:
    --------
    precursor = 'GATAGACAAGGTAGGAGAAAATGACTCGAACGAATTAGAGGTAGAGATAGATATCTATTCTATATTGAGAAGAGATAGAATAGAATCTGTAAAACGAGAAAATAATAAAACGTTTAGAAAGAGATGGGGTCTTACAAGGTCAAGAAAAGGCCTTACAAGGTCAAGAAAA'
    seq2test_effseq = 'AAAACGTTTAGAAAGAGATGG'
    df = create_structure_features(precursor, seq2test_effseq)
    '''
    prefold, premfe, Best_miR, BestScore, anyhit_data, twohit_data = common

    # Effective sequence setup
    eff_seq = seq2test_effseq
    df = pd.DataFrame({'eff_seq': [eff_seq]})
    df = cmf.mers123(df)  # Generate sequence features

    data = {}
    data['premfe'] = premfe

    # Location and distance calculations
    start = dist_5p = precursor.find(eff_seq)
    stop = dist_5p + len(eff_seq)
    dist_3p = len(precursor) - int(dist_5p) - len(eff_seq)
    data['dist_5p'], data['dist_3p'] = dist_5p, dist_3p

    # Structure substring analysis
    a = max(start - 4, 0)
    b = min(stop + 4, len(precursor))
    substring_fold = prefold[a: b]
    ref_seq = precursor[a: b]

    # Structural metrics calculations
    data['paired_percentage'] = cmf.get_paired_percentage(substring_fold)
    data['length_longest_bulge'] = cmf.length_largest_bulge(substring_fold)
    data['length_longest_loop'] = cmf.length_longest_loop(substring_fold)
    data['longest_paired_length'] = cmf.length_largest_bracket_sequence(substring_fold)
    data['mircheck_conclu'], data['fback_start'], data['fback_stop'] = hp.call_mirCheck(prefold, start, stop)

    # Rolling average calculations
    for window in [3, 5, 7]:
        data['paired_roll' + str(window)] = cmf.get_paired_rolling_average(substring_fold, window)

    # Motif analysis
    for motif in ['AAA', 'TTT', 'CCC', 'GGG']:
        data['NBtriplet' + motif[0]] = cmf.number_of_motif(ref_seq, motif)

    # Further structural feature calculations
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

    # Make output dataframe
    new_df = pd.DataFrame([data])
    df = df.drop(['eff_seq'], axis=1)
    df = pd.concat([df, new_df], axis=1)
    return df

def get_siRNA_structure(name, precursor, DicerCall=21, tmpdir='../tmp/', species='ath', mirbase_file='../dbs/mature.fa'): 
  if not os.path.exists(tmpdir): os.makedirs(tmpdir)
  timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

  prefold, premfe = ret.run_RNAfold(precursor)
  Best_miR, BestScore, anyhit_data, twohit_data = miRNA_target_search(precursor, tmpdir, species, mirbase_file)
  common = [prefold, premfe, Best_miR, BestScore, anyhit_data, twohit_data]

  df = pd.DataFrame()
  for i in range(len(precursor)-DicerCall+1):
    for j in range(DicerCall-2, DicerCall+3): #considering drift, ie. [19, 20, 21, 22, 23]
      sseq = precursor[i:i+j]
      df_tmp = create_structure_features(precursor, sseq, common)
      df_tmp['_longeur_'] = j
      df = pd.concat([df, df_tmp])
  siRNA_structure_file = tmpdir + timestamp + name + '_siRNA_structures.tsv'
  df.to_csv(siRNA_structure_file, sep='\t', index=False)
  return siRNA_structure_file

def run_one_precursor(name, precursor, DicerCall=21, 
                        tmpdir='../tmp/', species='ath', mirbase_file='../dbs/mature.fa'):
  siRNA_structure_file = get_siRNA_structure(name, precursor, DicerCall, tmpdir, species, mirbase_file)
  position_list = mlloc.classify_a_file(siRNA_structure_file)
  te, NB_instances = pretreat_location_features(siRNA_structure_file)
  siRNA_indication_file = encode_and_compute_weight(te, NB_instances, datafile)
  start, end, score, outfile, outfile2 = compute_indications_for_effector_start_end(siRNA_indication_file, position_list)
  print(f"{name}, predicted siRNA start: {start}, end: {end}, score: {score}")
  #bpi.draw_6candidates_interface(outfile, outfile2, precursor) #241015 tmporay add this line
  return start, end, score, outfile, outfile2

def dna_to_rna(dna_sequence):
    rna_sequence = dna_sequence.replace('T', 'U').replace('t', 'u')   
    return rna_sequence

def rna_to_dna(rna_sequence):
    dna_sequence = rna_sequence.replace('U', 'T').replace('u', 't')   
    return dna_sequence

def user_interface(name, pri, DicerCall, outdir):
    start, end, score, outfile, outfile2 = run_one_precursor(name, pri, DicerCall)
    outfile_newname = outdir + name + '.effector_localization_indication.tsv'
    outfile2_newname = outdir + name + '.effector_localization_top6_recommendation.tsv'
    os.rename(outfile, outfile_newname)
    os.rename(outfile2, outfile2_newname)
    bpi.draw_6candidates_interface(outfile_newname, outfile2_newname, pri)

if __name__ == "__main__":
  timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
  outdir = '../output/' + timestamp + '/'; os.makedirs(outdir)
  precursorName='yourPrecursor'
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
  pri = rna_to_dna(pri) # in case user provides RNA, convert it to DNA
  user_interface(precursorName, pri, DicerCall, outdir)
  print('==== See results in', outdir)
  pass
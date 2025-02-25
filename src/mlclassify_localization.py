'''
Chao-Jung Wu
2024-May-03
Update: 240729
update: 241008
Update: 241009

pickle
classify a file for siRNA localization
'''
import warnings
warnings.filterwarnings('ignore')
import random, sys
import pandas as pd
import numpy as np
from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder, QuantileTransformer
import pickle

def set_seed(seed=0):
  random.seed(seed)
  np.random.seed(seed)
seed = 0
set_seed(seed=seed)


file_used_as_background = '../dbs/background.tsv'
datafile = '../model/Arabidopsis_strcture_feature_importance_n_correlation.tsv' 
pickle_file='../model/RFAs100.pkl'#; print('Using RFAs100 model')


k = 100
tag='consistent'
cols_to_drop = ['CONTIG', 'eff_seq', 'retained', tag, 'segment']
cols2drop_for_locapredi = ['CONTIG', 'k', 'n', 'N', 'eff_strand', 'eff_pos', 'eff_frq', 'ext_k', 'length', 
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
show_cols = 'dist_5p, dist_3p, _longeur_, Probability_FALSE, Probability_TRUE, Predicted_Class'.split(', ')



def concatenate(df1, df2):
    result_df = pd.concat([df1, df2], ignore_index=True)
    result_df = result_df.reindex(columns=df2.columns)
    result_df = result_df.fillna(0)
    return result_df

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
 
def get_top_k_features(chi2_scores_dict, k=100):
    sorted_features = sorted(chi2_scores_dict.items(), key=lambda item: item[1], reverse=True)
    top_k_features = dict(sorted_features[:k])
    return top_k_features

def get_data(datafile, k=100):
   df = pd.read_csv(datafile, sep='\t')
   feature_importances = dict(zip(df['feature'], df['importance(chi2)']))
   feature_correlations = dict(zip(df['feature'], df['correlation(pointbiserialr)']))
   feature_names_topk = get_top_k_features(feature_importances, k).keys()
   return feature_names_topk, feature_correlations, feature_importances

def pretreat_location_features(infile_of_Instance):
    outfile = infile_of_Instance[:-3] + 'pretreated.tsv'

    df1 = pd.read_csv(infile_of_Instance, sep = '\t', index_col=False)
    NB_instances = len(df1)
    df1 = df1.drop(columns=cols_to_drop,axis=1, errors='ignore')
    df1['CONTIG'] = 'nana'
    df1['retained'] = True
    df1['consistent'] = True  

    df2 = pd.read_csv(file_used_as_background, sep = '\t', index_col=False)
    df2 = df2.drop(columns=cols_to_drop,axis=1, errors='ignore')
    df2['_longeur_'] = 0
    result_df = concatenate(df1, df2)
    result_df.to_csv(outfile, sep='\t', index=False)
    return outfile, NB_instances

def classify_a_file(infile2, pickle_file=pickle_file, datafile=datafile, cols_to_drop=cols_to_drop):
  infile2, NB_instances = pretreat_location_features(infile2)
  feature_names_topk, feature_correlations, feature_importances = get_data(datafile, k)  
  output_file = infile2[:-3] + 'prediction.tsv'
  
  # Load the instances to be predicted
  df = pd.read_csv(infile2, sep = '\t', index_col=False); #print(df.columns.values.tolist()) # Show list of features
  X = df.drop(columns=cols_to_drop,axis=1, errors='ignore')
  ord_col = []
  cat_col = [i for i, var in enumerate(X.columns) if X[var].dtype == 'O' or X[var].dtype == 'bool']
  num_col = [i for i, var in enumerate(X.columns) if X[var].dtype != 'O' and X[var].dtype != 'bool']
  X_encoded = preprocess_direct(X, num_col, cat_col, ord_col)
  X100 = X_encoded[feature_names_topk]
  print('NB of attributes after encoding:', len(X_encoded.columns))
  print('NB of attributes used for prediction:', len(X100.columns))

  # Remove background
  X100 = X100.head(NB_instances); #print(len(df))
  df = df.head(NB_instances)

  with open(pickle_file, 'rb') as f:
    pickle_model = pickle.load(f)
  predicted_classes = pickle_model.predict(X100)
  probabilities = pickle_model.predict_proba(X100)
  df['Probability_FALSE'] = probabilities[:, 0]
  df['Probability_TRUE'] = probabilities[:, 1]
  df['Predicted_Class'] = predicted_classes
  df = df[show_cols]
  df.to_csv(output_file, sep='\t', index=False)
  print("Output prediction file: \n", output_file)

  df['position'] = df['dist_5p'] + 1
  position_list = df[df['Predicted_Class'] == True]['position'].tolist()
  return list(set(position_list))






if __name__ == "__main__":
  if len(sys.argv) < 2:
    print("usage: python siWalk_classify_precursors.py $annotation_file $ref $pickle_file")
    print("Note: Ensure that ref and pickle_file are consistent, as they must be made from the same run.")
    print("example usage:")
    print("annotation_file=../UnitTest/20241010_022147yourPrecursor_siRNA_structures.tsv")
    print("ref=../model/test_feature_importance_n_correlation.tsv")
    print("pickle_file=../model/test_model.pkl")
    print("python siWalk_classify_localization.py $annotation_file")
    print("python siWalk_classify_localization.py $annotation_file $ref $pickle_file\n")
    print("Existing...")
    sys.exit(0)
  siRNA_structure_file=sys.argv[1]
  mydatafile = sys.argv[2] if len(sys.argv) > 2 else datafile
  mypickle_file = sys.argv[3] if len(sys.argv) > 3 else pickle_file



  print("Predicting siRNA localization from annotation file: \n", siRNA_structure_file)
  position_list = classify_a_file(siRNA_structure_file)
  print(position_list)
  pass
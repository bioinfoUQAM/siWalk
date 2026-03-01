'''
Train a GradientBoosting+AdaBoost model using all feature types (Module X).

Reads an annotated training TSV, encodes features, selects the top-100
attributes by chi2 importance, trains a GradientBoostingWithAdaBoost
classifier (default) or RandomForestWithAdaBoost, and writes the pickled
model and feature evaluation file to output/.
The resulting model (renamed GBA100.pkl) is intended for precursor prediction.

Author: Chao-Jung Wu
Date:   2024-May-03
'''
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import random, sys, os
from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder, QuantileTransformer
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pointbiserialr
from datetime import datetime
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

outdir = '../output/'
if not os.path.exists(outdir): os.makedirs(outdir)

def set_seed(seed=0):
  """Set Python and NumPy random seeds for reproducibility."""
  random.seed(seed)
  np.random.seed(seed)
seed = 0
set_seed(seed=seed)

tag = 'consistent'
k   = 100
cols_to_drop = ['CONTIG', 'eff_seq', 'retained', tag, 'segment']

def preprocess_direct(X, num_col, cat_col, ord_col):
  """Encode and scale features: quantile-transform numerical, one-hot categorical."""
  if all(isinstance(col, int) for col in num_col): num_col = X.columns[num_col]
  if all(isinstance(col, int) for col in cat_col): cat_col = X.columns[cat_col]
  if all(isinstance(col, int) for col in ord_col): ord_col = X.columns[ord_col]

  scaler     = QuantileTransformer(n_quantiles=10, random_state=seed)
  cat_encode = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
  ord_encode = OrdinalEncoder()

  X[num_col] = scaler.fit_transform(X[num_col])
  encoded_cat_array = cat_encode.fit_transform(X[cat_col])
  encoded_cat_df = pd.DataFrame(encoded_cat_array, columns=cat_encode.get_feature_names_out(cat_col), index=X.index)
  X = X.drop(columns=cat_col)
  X = pd.concat([X, encoded_cat_df], axis=1)
  X[ord_col] = ord_encode.fit_transform(X[ord_col])
  return X

def calculate_feature_importances_with_Gini_RF(X, y, model=RandomForestClassifier()):
  """Calculate feature importances using Random Forest Gini impurity."""
  model.fit(X, y)
  importances = model.feature_importances_
  return dict(zip(X.columns, importances))

def calculate_feature_importances_with_chi2(X_encoded, y):
  """Calculate feature importances using chi2 test."""
  chi2_selector = SelectKBest(chi2, k='all')
  chi2_selector.fit(X_encoded, y)
  return dict(zip(chi2_selector.feature_names_in_, chi2_selector.scores_))

def get_top_k_features(chi2_scores_dict, k=100):
  """Return the top-k features sorted by score (descending)."""
  sorted_features = sorted(chi2_scores_dict.items(), key=lambda item: item[1], reverse=True)
  return dict(sorted_features[:k])

def calculate_feature_correlations(X, class_labels):
  """
  Calculate point biserial correlation of each feature with the class label.

  Range: -1 to 1. Used when the class variable is binary.
  """
  feature_correlations = {}
  for feature in X.columns:
    correlation, _ = pointbiserialr(X[feature], class_labels)
    feature_correlations[feature] = correlation
  return feature_correlations

def print_importance_and_correlation_as_table_and_figure(filename, feature_importances, feature_correlations, method='gini'):
  """Write feature importance and correlation to TSV and produce a bar chart."""
  outfile = filename[:-3] + timestamp + '_feature_importance_n_correlation.tsv'
  with open(outfile, 'w') as fho:
    print('\t'.join(['feature', 'importance(' + method + ')', 'correlation(pointbiserialr)']), file=fho)
    for k, v in feature_importances.items():
      print('\t'.join([str(x) for x in [k, v, feature_correlations[k]]]), file=fho)
  vis_Feature_Importance_and_Correlation(outfile, method)
  return outfile

def vis_Feature_Importance_and_Correlation(infile, method='gini'):
  """Plot feature importance and correlation as a dual-axis horizontal bar chart."""
  col_importance  = 'importance(' + method + ')'
  col_correlation = 'correlation(pointbiserialr)'
  df = pd.read_csv(infile, sep='\t', index_col=False).sort_values(by=col_importance, ascending=True)
  y_pos = np.arange(len(df['feature']))
  ranks = np.arange(1, len(df['feature']) + 1)
  rank_labels = [str(rank) if rank == 1 or rank % 100 == 0 else '' for rank in ranks]

  fig, ax1 = plt.subplots(figsize=(10, 6))
  ax1.barh(y_pos, df[col_importance], color='blue', label='Importance', height=1, alpha=0.3)
  ax1.set_xlabel('Importance (' + method + ')', fontsize=18)
  ax1.set_yticks(y_pos)
  ax1.set_yticklabels(rank_labels[::-1])
  ax1.tick_params(left=False)
  ax1.set_ylabel('Rank by Importance', fontsize=18)

  ax2 = ax1.twiny()
  ax2.barh(y_pos - 0.2, df[col_correlation], color='red', label='Correlation', height=1)
  ax2.set_xlabel('Correlation (Point Biserial)', fontsize=18)

  handles1, labels1 = ax1.get_legend_handles_labels()
  handles2, labels2 = ax2.get_legend_handles_labels()
  ax1.legend(handles=handles2 + handles1, labels=labels2 + labels1, loc='lower right', fontsize=14)
  plt.tight_layout()
  plt.savefig(infile[:-3] + 'barplot.png')

def get_feature_importances(X_encoded, y, clf, method):
  if method == 'chi2': return calculate_feature_importances_with_chi2(X_encoded, y)
  if method == 'gini': return calculate_feature_importances_with_Gini_RF(X_encoded, y, clf)

def train_on_a_file(infile, classifier_name='GradientBoostingWithAdaBoost', cols_to_drop=cols_to_drop, k=k, tag=tag):
  """Train the model on *infile* and return (pickle_path, feature_evaluation_path)."""
  filename = outdir + timestamp + "_" + os.path.basename(infile)

  if classifier_name == 'GradientBoostingWithAdaBoost':
    grabo = GradientBoostingClassifier(learning_rate=0.2, max_features='log2', min_samples_leaf=4, n_estimators=200, random_state=seed)
    clf   = AdaBoostClassifier(estimator=grabo, random_state=seed)
  if classifier_name == 'RandomForestWithAdaBoost':
    rf  = RandomForestClassifier(max_depth=30, max_features='log2', min_samples_leaf=6, n_estimators=200, random_state=seed)
    clf = AdaBoostClassifier(estimator=rf, random_state=seed)

  df = pd.read_csv(infile, sep='\t', index_col=False)
  X  = df.drop(columns=cols_to_drop, axis=1)
  y  = df[tag]
  print('NB of instances in training:', len(df))
  print(f"Number of occurrences where {tag} is True:  {df[tag].sum()}")
  print(f"Number of occurrences where {tag} is False: {(~df[tag]).sum()}")

  ord_col = []
  cat_col = [i for i, var in enumerate(X.columns) if X[var].dtype == 'O' or X[var].dtype == 'bool']
  num_col = [i for i, var in enumerate(X.columns) if X[var].dtype != 'O' and X[var].dtype != 'bool']
  X_encoded = preprocess_direct(X, num_col, cat_col, ord_col)
  NBattribute = len(X_encoded.columns)
  print('NB of features after encoding:', NBattribute)

  feature_correlations = calculate_feature_correlations(X_encoded, y)
  method = 'chi2'  # method: {chi2, gini}; gini requires a forest-based classifier
  feature_importances = get_feature_importances(X_encoded, y, clf, method=method)
  datafile = print_importance_and_correlation_as_table_and_figure(filename, feature_importances, feature_correlations, method=method)

  if k == 'all': k = min(NBattribute * 0.8, 1000)
  feature_names_topk = get_top_k_features(feature_importances, k=k).keys()
  X100 = X_encoded[feature_names_topk]
  clf.fit(X100, y)

  pickle_file = outdir + timestamp + '_' + classifier_name + '_model.pkl'
  with open(pickle_file, 'wb') as f:
    pickle.dump(clf, f)
  return pickle_file, datafile


if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("example usage:")
    print("annotated_training_file=../dbs/background.tsv")
    print("python siWalk_pickle_precursor.py $annotated_training_file [$classifier]")
    print("# Optional classifier: {GradientBoostingWithAdaBoost, RandomForestWithAdaBoost}")
    print("# Default: GradientBoostingWithAdaBoost")
    sys.exit(0)
  infile = sys.argv[1]
  classifier_name = sys.argv[2] if len(sys.argv) > 2 else 'GradientBoostingWithAdaBoost'
  print('Timestamp:', timestamp)
  print('Classifier name:', classifier_name)
  pickle_file, datafile = train_on_a_file(infile, classifier_name)
  print('Trained model:', pickle_file)
  print('Feature evaluation file:', datafile)
  pass

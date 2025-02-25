'''
Chao-Jung Wu
2024-May-03
Update: 240729
Update: 241009
Update: 241106: default build RFAs100 model
pickle precursor
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

outdir='../output/'
if not os.path.exists(outdir): os.makedirs(outdir)

def set_seed(seed=0):
  ''' random comes from many places. Here my pipeline could be affected by python and array randomness. Set seed to make the experiments reproducible. '''
  random.seed(seed) #python's randomness
  np.random.seed(seed) # random in arrays
seed = 0
set_seed(seed=seed)

tag='consistent'
k=100
cols_to_drop = ['CONTIG', 'eff_seq', 'retained', tag, 'segment'] #, 'group']

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
 
def calculate_feature_importances_with_Gini_RF(X, y, model=RandomForestClassifier()):
    """
    Calculate feature importances using RandomForestClassifier.
    Gini Importance (Mean Decrease in Impurity)
    Gini impurity is a measure of how often a randomly chosen element from the dataset would be incorrectly labeled if it were randomly labeled according to the distribution of labels in the subset.
    Args:
        X (pd.DataFrame): The feature matrix.
        y (pd.Series): The target vector (class labels).
    Returns:
        feature_importances (dict): A dictionary where keys are feature names and values are their importances. """
    model.fit(X, y)
    importances = model.feature_importances_
    feature_importances = dict(zip(X.columns, importances))
    return feature_importances

def calculate_feature_importances_with_chi2(X_encoded, y):
    chi2_selector = SelectKBest(chi2, k='all')
    chi2_selector.fit(X_encoded, y)
    chi2_scores = chi2_selector.scores_
    feature_names = chi2_selector.feature_names_in_
    chi2_scores_dict = dict(zip(feature_names, chi2_scores))
    return chi2_scores_dict

def get_top_k_features(chi2_scores_dict, k=100):
    sorted_features = sorted(chi2_scores_dict.items(), key=lambda item: item[1], reverse=True)
    top_k_features = dict(sorted_features[:k])
    return top_k_features

def calculate_feature_correlations(X, class_labels):
    """
    Calculate the point biserial correlation of each feature with the class column.
    Measures correlation between a binary and a continuous variable.
    Range: -1 to 1.
    Usage: Specifically used when one variable is binary.
    Args:
        X (pd.DataFrame): The dataframe containing features but not the class column.
        class_labels: The dataframe representing the class labels.
    Returns:
        feature_correlations (dict): A dictionary where keys are feature names and values are correlations.
    """
    feature_correlations = {}
    for feature in X.columns:
        correlation, _ = pointbiserialr(X[feature], class_labels)
        feature_correlations[feature] = correlation  
    return feature_correlations

def print_importance_and_correlation_as_table_and_figure(filename, feature_importances, feature_correlations, method='gini'):
  outfile = filename[:-3] + timestamp + '_feature_importance_n_correlation.tsv'
  fho = open (outfile, 'w')
  print('\t'.join(['feature', 'importance('+method+')', 'correlation(pointbiserialr)']), file=fho)
  for k, v in feature_importances.items():
      data = [k, v, feature_correlations[k]]
      print('\t'.join([str(x) for x in data]), file=fho)
  fho.close()
  vis_Feature_Importance_and_Correlation(outfile, method)
  return outfile

def vis_Feature_Importance_and_Correlation(infile, method='gini'):
    col_importance = 'importance('+method+')'
    col_correlation = 'correlation(pointbiserialr)'
    df = pd.read_csv(infile, sep='\t', index_col=False)
    df = df.sort_values(by=col_importance, ascending=True)
    
    # Generate y positions for the features (reversed for top-down ranking)
    y_pos = np.arange(len(df['feature']))
    ranks = np.arange(1, len(df['feature']) + 1)
    
    # Filter ranks to only show multiples of 100 (e.g., 1, 100, 200, 300, ...)
    rank_labels = [str(rank) if rank == 1 or rank % 100 == 0 else '' for rank in ranks]
    
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.barh(y_pos, df[col_importance], color='blue', label='Importance', height=1, alpha=0.3)
    ax1.set_xlabel('Importance ('+method+')', fontsize=18)   
    # ax1.set_xlim(0, 0.05)  # Adjust the upper limit as needed to better display the smaller values
    ax1.set_yticks(y_pos)  # Set y positions
    ax1.set_yticklabels(rank_labels[::-1])  # Reverse rank labels to show top-down
    ax1.tick_params(left=False)  # Hide the tick marks
    ax1.set_ylabel('Rank by Importance', fontsize=18)

    # Plot correlation scores on the same figure, but secondary axis
    ax2 = ax1.twiny()
    ax2.barh(y_pos - 0.2, df[col_correlation], color='red', label='Correlation', height=1)
    ax2.set_xlabel('Correlation (Point Biserial)', fontsize=18)
    # ax2.set_xlim(-.2, .2)
    
    # Combine the legends of both axes into one
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles=handles2 + handles1, labels=labels2 + labels1, loc='lower right', fontsize=14)

    plt.tight_layout()
    plt.savefig(infile[:-3] + 'barplot.png')

def get_feature_importances(X_encoded, y, clf, method):
  if method == 'chi2':
    feature_importances = calculate_feature_importances_with_chi2(X_encoded, y)
  if method == 'gini':
    feature_importances = calculate_feature_importances_with_Gini_RF(X_encoded, y, clf)
  return feature_importances

def train_on_a_file(infile, classifier_name='GradientBoostingWithAdaBoost', cols_to_drop=cols_to_drop, k=k, tag=tag):
  filename = outdir + timestamp + "_" + os.path.basename(infile)

  if classifier_name == 'GradientBoostingWithAdaBoost':
    grabo = GradientBoostingClassifier(learning_rate=0.2, max_features='log2', min_samples_leaf=4, n_estimators=200, random_state=seed)
    clf = grabo_ada = AdaBoostClassifier(estimator=grabo, random_state=seed)
  if classifier_name == 'RandomForestWithAdaBoost':
    rf = RandomForestClassifier(max_depth=30, max_features='log2', min_samples_leaf=6, n_estimators=200, random_state=seed)
    clf = rf_ada = AdaBoostClassifier(estimator=rf, random_state=seed)
  
  # Load the instances to be learned
  df = pd.read_csv(infile, sep = '\t', index_col=False)#; print(df.head())
  X = df.drop(columns=cols_to_drop,axis=1)
  y = df[tag]
  print('NB of instances in training:', len(df))
  count_true = df[tag].sum()
  count_false = (~df[tag]).sum()
  print(f"Number of occurrences where {tag} is True: {count_true}")
  print(f"Number of occurrences where {tag} is False: {count_false}")

  ord_col = []
  cat_col = [i for i, var in enumerate(X.columns) if X[var].dtype == 'O' or X[var].dtype == 'bool']
  num_col = [i for i, var in enumerate(X.columns) if X[var].dtype != 'O' and X[var].dtype != 'bool']
  X_encoded = preprocess_direct(X, num_col, cat_col, ord_col)
  NBattribute = len(X_encoded.columns)
  print('NB of features after encoding:', NBattribute)

  feature_correlations = calculate_feature_correlations(X_encoded, y)
  method = 'chi2' # method={chi2, gini}; with gini, clf has to be forest-based.
  feature_importances = get_feature_importances(X_encoded, y, clf, method=method)
  datafile = print_importance_and_correlation_as_table_and_figure(filename, feature_importances, feature_correlations, method=method)

  # Get the best k feature names
  if k=='all': k=min(NBattribute*0.8, 1000)
  top_k_features = get_top_k_features(feature_importances, k=k)
  feature_names_topk = top_k_features.keys()
  # Train using only subset best features to obtain the desired model
  X100 = X_encoded[feature_names_topk]
  clf.fit(X100, y)
  # Pickle the trained model
  pickle_file = outdir + timestamp + '_' + classifier_name + '_model.pkl'
  with open(pickle_file, 'wb') as f:
    pickle.dump(clf, f)
  return pickle_file, datafile




if __name__ == "__main__":
  if len(sys.argv) != 2:
      print("example usage:")
      print("annotated_training_file=../dbs/background.tsv")
      print("python siWalk_pickle_precursor.py $annotated_training_file [$classifier]")
      print("# Defining classifier is optional, options are {GradientBoostingWithAdaBoost, RandomForestWithAdaBoost}, default is GradientBoostingWithAdaBoost")
      sys.exit(0)
  infile = sys.argv[1]
  classifier_name = sys.argv[2] if len(sys.argv) > 2 else 'GradientBoostingWithAdaBoost'
  print('Timestamp:', timestamp)
  print('Classifier name:', classifier_name)
  pickle_file, datafile = train_on_a_file(infile, classifier_name)
  print('Trained model:', pickle_file)
  print('Feature evaluation file:', datafile)
  pass


  # infile = '../input/merged12_file.tsv'
  # infile = '../input/training_set_all29libs.tsv'
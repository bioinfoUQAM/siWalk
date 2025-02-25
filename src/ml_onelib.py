'''
Chao-Jung Wu
2024-May-03
Update: 2024-July-26

ML pipeline
feature selection: chi2
model: random forest, GradientBoosting, LogisticRegression, svm
select best 100 features with higher chi2

Take one library's contig features which non-features are removed.


## to kill a job
ps -ef
kill <PID>
'''
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import random

from sklearn.model_selection import train_test_split, cross_val_score, cross_validate
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OrdinalEncoder, MinMaxScaler, OneHotEncoder, QuantileTransformer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, chi2, mutual_info_classif
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score, precision_score
from sklearn.metrics import recall_score, balanced_accuracy_score, roc_auc_score, r2_score
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.ensemble import GradientBoostingClassifier

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import plot_bar_with_std



def set_seed(seed=0):
  ''' random comes from many places. Here my pipeline could be affected by python and array randomness. Set seed to make the experiments reproducible. '''
  random.seed(seed) #python's randomness
  np.random.seed(seed) # random in arrays

# Scoring methods
from sklearn.metrics import matthews_corrcoef, make_scorer
mcc_scorer = make_scorer(matthews_corrcoef)
def tp_scorer(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    return cm[1, 1]  # True Positives
def tn_scorer(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    return cm[0, 0]  # True Negatives
def fp_scorer(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    return cm[0, 1]  # False Positives
def fn_scorer(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    return cm[1, 0]  # False Negatives
def specificity_(y_true, y_pred):
    tn = np.sum((y_true == 0) & (y_pred == 0))
    fp = np.sum((y_true == 0) & (y_pred == 1))
    return tn / (tn + fp)
tp_scorer = make_scorer(tp_scorer)
tn_scorer = make_scorer(tn_scorer)
fp_scorer = make_scorer(fp_scorer)
fn_scorer = make_scorer(fn_scorer)
specificity_scorer = make_scorer(specificity_)
''' scoring is a list if using only built-in scorer #https://scikit-learn.org/stable/modules/model_evaluation.html
    if incoporating make_scorer, need to set up a dictionary
scoring = ['precision', 'recall', 'f1_weighted', 'f1', 'roc_auc', 'accuracy'] '''
scoring = {
    'tp': tp_scorer,
    'tn': tn_scorer,
    'fp': fp_scorer,
    'fn': fn_scorer,
    'precision': 'precision',
    'sensitivity': 'recall',
    'specificity': specificity_scorer,
    'acc': 'accuracy',
    'mcc': mcc_scorer, 
    'auc': 'roc_auc',
    'f1': 'f1',
}



seed = 0
set_seed(seed=seed)
cv = 5 #cross validation folds
refit = 'precision' # score to optimize

# Feature selection meethods
fs_methods = {'chi2': chi2}
k_best_features = [100]

# Held One Out CV Classifiers (parameters chosen by grid search)
rf = RandomForestClassifier(max_depth=10, max_features='log2', min_samples_leaf=4, min_samples_split=10, n_estimators=300, random_state=seed)
grabo = GradientBoostingClassifier(criterion='squared_error', max_depth=7, max_features='log2', random_state=seed, subsample=0.8)
lr = LogisticRegression(C=0.1, max_iter=10000, random_state=seed, solver='sag')
svcrbf = svm.SVC(C=1, degree=2, gamma='auto', random_state=seed, probability=True)
clfs = classifiers = [grabo, rf, lr, svcrbf]

def Matthews_correlation_coefficient (tp, tn, fp, fn):
  ''' same as in MirDup, Leclercq, NAR 2013
  When the positive and negative classes are unbalanced, F1 can not be used. Instead, use MCC.
  The MCC value is between -1 and 1.
  0 : random
  -1: disagree
  1 : agree '''
  from math import sqrt
  up = tp*tn - fp*fn
  dn = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if dn == 0: return 0
  return up/dn

def specificity (tn, fp):
  ''' same as in MirDup, Leclercq, NAR 2013 '''
  if tn + fp == 0: return 0
  return tn/(tn + fp)

def evaluate_prediction_performance (y_true, y_pred):
  ''' from sklearn.metrics import accuracy_score, f1_score, recall_score, r2_score, etc
  https://scikit-learn.org/stable/modules/model_evaluation.html
  @balanced_accuracy_score: accuracy of imbalanced datasets. Average of recall obtained on each class.
  @f1 weighted: Account for label imbalance
  @f1 binary: Applicable only if targets are binary
  @f1 samples: For multilabel classification
  @f1 macro: Does not take label imbalance into account; irrelevant to my problem
  '''
  #print('========== Prediction performance =============')
  cm = confusion_matrix(y_true, y_pred)
  [tn, fp], [fn, tp] = cm[0], cm[1]

  #= literature metrics
  mcc = Matthews_correlation_coefficient (tp, tn, fp, fn)
  speci = specificity (tn, fp)
  
  #= sklearn.metrics
  acc = accuracy_score(y_true, y_pred)
  f1b = f1_score(y_true, y_pred, average='binary')
  preci = precision_score(y_true, y_pred)
  sens = recall_score(y_true, y_pred) #aka sensitivity or true positive rate
  auc = roc_auc_score(y_true, y_pred)
  
  title = 'tp, tn, fp, fn, precision, sensitivity, specificity, acc, mcc, auc, f1'.split(', ')
  print('\t'.join(title))
  data = [str(round(x,3)) for x in [tp, tn, fp, fn, preci, sens, speci, acc, mcc, auc, f1b]]
  line = '\t'.join(data);print(line)

def preprocessing (num_col, cat_col, ord_col):
  ''' Found unknown categories ['blah blah'] in column ### during transform 
  Solutions: (1) handle_unknown='ignore'; or (2) handle_unknown='use_encoded_value', unknown_value=-1 '''
  # Transformers for numerical data
  numerical_transformer = Pipeline(steps=[
        ('scaler', QuantileTransformer(n_quantiles=10, random_state=seed)),
  ])
  # Transformers for categorical data
  categorical_transformer = Pipeline(steps=[
        ("encoder", OneHotEncoder(categories='auto', handle_unknown='ignore')),
  ])
  # Transformers for ordinal data
  ordinal_transformer = Pipeline(steps=[
        ("encoder", OrdinalEncoder()),#(handle_unknown='use_encoded_value', unknown_value=-1)),
  ])

  # Combine pipelines using ColumnTransformer
  preproc_pipe = ColumnTransformer(
      transformers=[
        ('num', numerical_transformer, num_col),
        ('cat', categorical_transformer, cat_col),
        ('ord', ordinal_transformer, ord_col),
    ],
    remainder="passthrough"
  )
  return preproc_pipe

def build_complete_pipeline (num_col, cat_col, ord_col, fs, clf):
  preproc_pipe = preprocessing (num_col, cat_col, ord_col)
  complete_pipe = Pipeline(steps=[
        ("preprocessor", preproc_pipe),
        ("feature_selection", fs),
        ("classifier", clf),
  ])
  return complete_pipe

def build_short_pipeline(fs, clf):
  short_pipe = Pipeline(steps=[
        ("feature_selection", fs),
        ("classifier", clf),
  ])
  return short_pipe

def preprocess_direct(X, num_col, cat_col, ord_col):
  scaler = QuantileTransformer(n_quantiles=10, random_state=seed)
  cat_encode = OneHotEncoder(handle_unknown='ignore')
  ord_encode = OrdinalEncoder()
  X.iloc[:,num_col] = scaler.fit_transform(X.iloc[:,num_col])
  X.iloc[:,cat_col] = cat_encode.fit_transform(X.iloc[:,cat_col]).shape[0]
  X.iloc[:,ord_col] = ord_encode.fit_transform(X.iloc[:,ord_col])
  return X

def display_feature_selection (X, y, fs, name, infile):
  ''' eg: fs = SelectKBest(score_func=mutual_info_classif) '''
  outfile = infile[:-3] + 'report_HeldOneOut.feature_rank_scores.' + name + '.tsv'
  fho = open (outfile, 'w')
    
  outfig = infile[:-3] + name + '_feature_ranking.png'
  fs.fit(X, y)
  features = fs.feature_names_in_
  scores = fs.scores_
  d = dict(sorted(zip(features, scores), key=lambda item: item[1]))

  print('#====' + name + ' feature selection', file=fho)
  for k, v in d.items(): print(k, 'feature score:\t', round(v, 10), file=fho)
  fho.close()

  import matplotlib.pyplot as plt
  n = round(len(d)/30) + 1
  plt.figure(figsize=(8, 6*n)) #figsize=(wide, long)

  plt.title(name + ' feature selection')
  plt.xlabel('score (the higher the better)')
  plt.yticks(rotation=45)
  plt.barh(list(d.keys()), list(d.values()))
  plt.savefig(outfig)
  print('\nsee feature ranking png:', outfig)
  print('see feature_rank_scores:', outfile)

def demo_train_test_split_once(X_train,y_train, X_test, y_test, k_best_features, num_col, cat_col, ord_col, infile):
  print('Split train and test sets, no cross validation')
  outfig = infile[:-3] + 'roc_curve.png'
  plt.figure()
  color = cm.tab20(np.linspace(0, 1, 50))
  count = 0
   
  for clf in clfs:
    print('\n\n===== classifier:', clf)
    for name, method in fs_methods.items():
      for k in k_best_features:
        fs = SelectKBest(method, k=k)
        print('========= feature selection:', fs, '\n')
        complete_pipe = build_complete_pipeline (num_col, cat_col, ord_col, fs, clf)
        complete_pipe #Visualize pipeline in jupytre
        # Run ML pipeline
        results = complete_pipe.fit(X_train,y_train)
        # Evaluate pipeline accuracy etc
        predictions = complete_pipe.predict(X_test)
        evaluate_prediction_performance (y_test, predictions)

        y_score = complete_pipe.predict_proba(X_test)[:, 1]
        # Compute ROC curve and ROC area
        fpr, tpr, _ = roc_curve(y_test, y_score)
        roc_auc = auc(fpr, tpr)
        c = color[count]; count += 1
        plt.plot(fpr, tpr, color=c, lw=2, label=clf.__class__.__name__ + ' (area = %0.2f)' % roc_auc)
  plt.xlabel('False Positive Rate')
  plt.ylabel('True Positive Rate')
  plt.title('\nReceiver Operating Characteristic')
  plt.legend(loc="lower right")
  plt.savefig(outfig)
  print('see Receiver Operating Characteristic png:', outfig)

def cross_validation_Held_One_Out__with_display_feature_selection (X, y, k_best_features, num_col, cat_col, ord_col, infile):
  outfile = infile[:-3] + 'report_HeldOneOut.tsv'
  fho = open (outfile, 'w')
 
  keys = scoring.keys()
  xxx, yyy = ', '.join(keys), '_std, '.join(keys) + '_std'
  title = '#model, model_combination, , '+ xxx + ', , ' + yyy
  print('\t'.join(title.split(', ')), file=fho)

  print('Cross validation strategy = Held_One_Out;\t cross validation folds = ', cv)
  count = 1
  for fsname, method in fs_methods.items():
    for clf in clfs:
      print('\n\n===== classifier:', clf)
      for k in k_best_features:
        fs = SelectKBest(method, k=k)
        print('========= feature selection:', fs)
        # Run ML pipeline
        X = preprocess_direct(X, num_col, cat_col, ord_col)
        short_pipe = build_short_pipeline(fs, clf)
        short_pipe #Visualize pipeline in jupytre
        cv_result = cross_validate(short_pipe, X, y, cv=5, scoring=scoring)
        print(cv_result)
        # Evaluate pipeline accuracy etc
        mylist_ave = ['average'] + [round(sum(cv_result['test_' + s]) / len(cv_result['test_' + s]), 3) for s in scoring]
        std_ = ['std'] + [round(np.std(cv_result['test_' + s]), 3) for s in scoring]
        model_combination = str(clf).replace("\n", "").replace(" ", "") + '.' + fsname + '.k' + str(k)
        mylist_ave_std = [count, model_combination] + mylist_ave + std_
        print('\t'.join([str(x) for x in mylist_ave_std]), file=fho, flush=True)
        count += 1
    display_feature_selection (X, y, fs, fsname, infile)
  fho.close()
  outfig = plot_bar_with_std.run(outfile)
  print('\n\nsee report file for performance evaluation:', outfile)
  print('see report png:', outfig)
  return outfile, outfig

def prepare_dataset(infile, k_best_features, classtag='retained'):
    df = pd.read_csv(infile, sep = '\t', index_col=False)
    # X features; y class labels 
    X = df.drop([classtag],axis=1)
    y = df[classtag]
    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=seed)
    # Identify numerical and categorical columns
    ord_col = []
    cat_col = [i for i, var in enumerate(X.columns) if X[var].dtype == 'O' or X[var].dtype == 'bool']   
    num_col = [i for i, var in enumerate(X.columns) if X[var].dtype != 'O' and X[var].dtype != 'bool']
    # If k_best_features = ['all', 100, 1000] but number of features only 300, then remove 1000, return k_best_features = ['all', 100]
    k_best_features = [x for x in k_best_features if isinstance(x, str) or x <= len(X.columns)+10]
    return X, y, X_train,y_train, X_test, y_test, num_col, cat_col, ord_col, k_best_features

def run(infile, k_best_features=[100]):
    X, y, X_train, y_train, X_test, y_test, num_col, cat_col, ord_col, k_best_features = prepare_dataset(infile, k_best_features)
    outfile, outfig = cross_validation_Held_One_Out__with_display_feature_selection(X, y, k_best_features, num_col, cat_col, ord_col, infile)
    return outfile, outfig

def demo_once(infile, k_best_features=[100]):
    X, y, X_train, y_train, X_test, y_test, num_col, cat_col, ord_col, k_best_features = prepare_dataset(infile, k_best_features)
    demo_train_test_split_once(X_train,y_train, X_test, y_test, k_best_features, num_col, cat_col, ord_col, infile)

def test2(infile, k_best_features=[100]):
    X, y, X_train, y_train, X_test, y_test, num_col, cat_col, ord_col, k_best_features = prepare_dataset(infile, k_best_features)
    X = preprocess_direct(X, num_col, cat_col, ord_col)
    print(X)
    X.to_csv('test.tsv', sep='\t', index = False)










if __name__ == "__main__":
  infile = '../UnitTest_ml/GSM1087987_C0FCSb.contig_features.rm_nonfeat.tsv'
  # outfile, outfig = run(infile, k_best_features=[100])
  
  test2(infile)
  pass
'''
Author: Chao-Jung Wu
Date: 2017-11-07
'''
from math import sqrt

def specificity (tn, fp):
  ''' same as in MirDup, Leclercq, NAR 2013 '''
  if tn + fp == 0: return 0
  return round (tn/(tn + fp), 3)

def precision (tp, fp):
  if tp + fp == 0: return 0
  return round (tp / (tp + fp), 3)

def sensitivity (tp, fn):
  ''' same as in MirDup, Leclercq, NAR 2013
  aka Recall '''
  if tp + fn == 0: return 0
  return round (tp / (tp + fn), 3)

def F1_score (tp, tn, fp, fn):
  if 2 * tp + fp + fn == 0: return 0
  f1 = 2 * tp / (2 * tp + fp + fn)
  return round (f1, 3)

def accuracy (tp, tn, fp, fn):
  ''' same as in MirDup, Leclercq, NAR 2013 '''
  if tp + tn + fp + fn == 0: return 0
  acc = (tp + tn)/ float(tp + tn + fp + fn)
  return round (acc, 3)

def Matthews_correlation_coefficient (tp, tn, fp, fn):
  ''' same as in MirDup, Leclercq, NAR 2013
  When the positive and negative classes are unbalanced, F1 can not be used. Instead, use MCC.
  
  The MCC value is between -1 and 1.
  0 : random
  -1: disagree
  1 : agree '''
  up = tp*tn - fp*fn
  dn = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if dn == 0: return 0
  return round (up/dn, 3)
  
def report (tp, tn, fp, fn, msg=''):
  f1 = F1_score (tp, tn, fp, fn)
  acc = accuracy (tp, tn, fp, fn)
  sens = sensitivity (tp, fn)
  preci = precision (tp, fp)
  mcc = Matthews_correlation_coefficient (tp, tn, fp, fn)
  spec = specificity (tn, fp)

  title = (', tp, tn, fp, fn, precision, sensitivity, specificity, acc, mcc, f1').split(', ')
  print('\t'.join(title))
  data = [str(x) for x in [msg, tp, tn, fp, fn, preci, sens, spec, acc, mcc, f1]]
  line = '\t'.join(data)
  print(line)
  return f1, acc, mcc, preci, sens, spec

def example_usage ():
  tp = 22
  fp = 0
  tn = 1000-fp
  fn = 100-tp
  print('tp, tn, fp, fn =', tp, tn, fp, fn)
  report (tp, tn, fp, fn)

  negset = 1838879 - 141  
  tp = 29
  fp = 343
  tn = negset-fp
  fn = 141-tp
  print('tp, tn, fp, fn =', tp, tn, fp, fn)
  report (tp, tn, fp, fn)
  
  tp, tn, fp, fn = 26, 2300, 120, 2
  print('tp, tn, fp, fn =', tp, tn, fp, fn)
  report (tp, tn, fp, fn)

#example_usage ()
'''
Nucleotide content and Argonaute preference features.

For each segment in the summary file, queries the sorted BAM file via samtools
to retrieve DicerCall-nt reads and computes:
  - Overall A/C/G/T% and GC% (expression-weighted).
  - Per-position nucleotide frequencies at the 5-prime end, 3-prime end,
    and mid-region of the active siRNA (5 positions each, 66 features total).

Author: Chao-Jung Wu
Date:   2024-Apr-04
'''
import pandas as pd
import subprocess as sbp
import os

def getRevComp(seq):
  """Return the reverse complement of *seq* (DNA, uppercase)."""
  intab, outab = "ACGT", "TGCA"
  trantab = str.maketrans(intab, outab)
  seq = seq.upper()
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

def cal_acgt_percent(total, a, c, g, t):
  """Return (A%, C%, G%, T%, GC%) rounded to 2 decimal places."""
  if total == 0: return 0, 0, 0, 0, 0
  ap = round(100 * a / total, 2)
  cp = round(100 * c / total, 2)
  gp = round(100 * g / total, 2)
  tp = round(100 * t / total, 2)
  gc_content = round(gp + cp, 2)
  return ap, cp, gp, tp, gc_content

def argonaute_preference(samdata, DicerCall):
  '''
  Compute expression-weighted nucleotide frequencies at key siRNA positions.

  For DicerCall-nt reads, records nucleotide counts at the first 5 positions
  (5-prime), last 5 positions (3-prime), and middle 5 positions of each read.
  Returns 66 features (6 global + 15 positions × 4 nucleotides).

  @samdata   : list of split SAM record lines.
  @DicerCall : expected siRNA length (e.g. 21).
  '''
  d = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_5p1 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_5p2 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_5p3 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_5p4 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_5p5 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_3p1 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_3p2 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_3p3 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_3p4 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_3p5 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_md1 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_md2 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_md3 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_md4 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
  d_md5 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

  for i in samdata:
    FLAG, SEQ, lenSEQ = i[1], i[9], len(i[9])

    if lenSEQ == DicerCall:
      for j in range(lenSEQ): d[SEQ[j]] += 1

      if FLAG == '16': SEQ = getRevComp(SEQ)
      elif FLAG == '0': pass
      else: print('error, strand not determined, FLAG =', FLAG)

      n5p1, n5p2, n5p3, n5p4, n5p5 = SEQ[0], SEQ[1], SEQ[3], SEQ[4], SEQ[5]
      n3p1, n3p2, n3p3, n3p4, n3p5 = SEQ[-1], SEQ[-2], SEQ[-3], SEQ[-4], SEQ[-5]
      mid1, mid2, mid3, mid4, mid5  = SEQ[9], SEQ[10], SEQ[11], SEQ[12], SEQ[13]

      d_5p1[n5p1] += 1; d_5p2[n5p2] += 1; d_5p3[n5p3] += 1
      d_5p4[n5p4] += 1; d_5p5[n5p5] += 1
      d_3p1[n3p1] += 1; d_3p2[n3p2] += 1; d_3p3[n3p3] += 1
      d_3p4[n3p4] += 1; d_3p5[n3p5] += 1
      d_md1[mid1] += 1; d_md2[mid2] += 1; d_md3[mid3] += 1
      d_md4[mid4] += 1; d_md5[mid5] += 1

  total_frq_dicercall = sum(d.values())
  ap, cp, gp, tp, gc_content = cal_acgt_percent(total_frq_dicercall, d['A'], d['C'], d['G'], d['T'])
  plenty_features = [total_frq_dicercall, ap, cp, gp, tp, gc_content]

  many = (d_5p1, d_5p2, d_5p3, d_5p4, d_5p5,
          d_3p1, d_3p2, d_3p3, d_3p4, d_3p5,
          d_md1, d_md2, d_md3, d_md4, d_md5)
  for d in many:
    d = normalize(d)
    for k, v in sorted(d.items()):
      plenty_features.append(v)

  return [str(x) for x in plenty_features]

def normalize(d):
  """Return A/C/G/T percentages (0–100) for a nucleotide count dictionary."""
  total = sum(d.values())
  d2 = {}
  for k, v in sorted(d.items()):
    if total == 0: d2[k] = 0
    else: d2[k] = round(100 * v / total, 2)
  return d2

def get_samdata_from_file(infile):
  """Read a SAM file and return its records as a list of split lines."""
  with open(infile, 'r') as fh:
    samdata = [x.rstrip('\n').split('\t') for x in fh.readlines()]
  return samdata

def title_of_66_features():
  """Return the ordered column names for the 66 Argonaute-preference features."""
  t = '5p1, 5p2, 5p3, 5p4, 5p5, 3p1, 3p2, 3p3, 3p4, 3p5, md1, md2, md3, md4, md5'.split(', ')
  x = 'A, C, G, T'.split(', ')
  title = 'total_frq_DicerCall, A%, C%, G%, T%, GC%'.split(', ')
  for i in t:
    title += [i + j for j in x]
  return title

def print_argonaute(data, summaryFile):
  """Append Argonaute-preference columns to summaryFile and return the new path."""
  with open(summaryFile, 'r') as fh:
    first = fh.readline().rstrip('\n')

  outfile = summaryFile[:-3] + 'argonautestat.tsv'
  with open(outfile, 'w') as fho:
    title = title_of_66_features()
    print(first + '\t' + '\t'.join(title), file=fho)
    for i in data:
      print('\t'.join(i), file=fho)
  return outfile


class argonaute_class():

  def __init__(self, infile, DicerCall, bamfile):
    self.DicerCall = DicerCall
    self.bamfile = bamfile
    with open(infile, 'r') as fh:
      L = fh.readline().rstrip('\n').split('\t')
    self.ich = L.index('chr')
    self.ilb = L.index('L_bound')
    self.irb = L.index('R_bound')

  def caculate_acgt_content_from_slice(self, e):
    """Query BAM for reads in the segment interval and compute Argonaute features."""
    ch, lb, rb = e[self.ich], e[self.ilb], e[self.irb]
    interval = str(ch) + ':' + str(lb) + '-' + str(rb)
    cmd = ['samtools', 'view', self.bamfile, interval]
    stdout = sbp.check_output(cmd)
    samdata = [x.split() for x in stdout.decode().rstrip().split('\n')]
    if samdata == [[]] : return interval, []
    plenty_features = argonaute_preference(samdata, self.DicerCall)
    return e + plenty_features


if __name__ == "__main__":
  pass

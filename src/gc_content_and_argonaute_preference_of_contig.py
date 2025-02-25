'''
Date: 2024-Apr-04
DicerCall=21 for example
In a contig, keep those 21nt RNAs, and ignore other sizes.
Do statistics on these 21nt species, inclduing the overall A% and GC%; and the positional 5p1A% etc.
'''
import pandas as pd
import subprocess as sbp
import os

def getRevComp (seq):
  ''' updated: 2024-Apr-26 '''
  intab, outab = "ACGT", "TGCA"
  trantab = str.maketrans(intab, outab)
  seq = seq.upper()
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

def cal_acgt_percent (total, a, c, g, t):
  if total == 0: return 0, 0, 0, 0, 0
  ap = round(100* a / total, 2)
  cp = round(100* c / total, 2)
  gp = round(100* g / total, 2)
  tp = round(100* t / total, 2)
  gc_content = round(gp + cp, 2)
  return ap, cp, gp, tp, gc_content

def argonaute_preference (samdata, DicerCall):
  '''
  For any cise of DicerCall, record each of the first 5 nucleotide frequencies, the last 5 nt freq, and the mid 5 nt freq.
  Also return the gc% of the contig
  return 61 features
  '''
  d = {'A':0, 'T':0, 'C':0, 'G':0}

  d_5p1 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_5p2 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_5p3 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_5p4 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_5p5 = {'A':0, 'T':0, 'C':0, 'G':0}

  d_3p1 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_3p2 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_3p3 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_3p4 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_3p5 = {'A':0, 'T':0, 'C':0, 'G':0}

  d_md1 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_md2 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_md3 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_md4 = {'A':0, 'T':0, 'C':0, 'G':0}
  d_md5 = {'A':0, 'T':0, 'C':0, 'G':0}  

  for i in samdata:
    FLAG, SEQ, lenSEQ = i[1], i[9], len(i[9])

    if lenSEQ == DicerCall:
      for j in range(lenSEQ): d[SEQ[j]] += 1

      if FLAG == '16': SEQ = getRevComp (SEQ)
      elif FLAG == '0': pass
      else: print('error, strand not determined, FLAG =', FLAG)

      n5p1, n5p2, n5p3, n5p4, n5p5 = SEQ[0], SEQ[1], SEQ[3], SEQ[4], SEQ[5]
      n3p1, n3p2, n3p3, n3p4, n3p5 = SEQ[-1], SEQ[-2], SEQ[-3], SEQ[-4], SEQ[-5]
      mid1, mid2, mid3, mid4, mid5 = SEQ[9], SEQ[10], SEQ[11], SEQ[12], SEQ[13]

      d_5p1[n5p1] += 1
      d_5p2[n5p2] += 1
      d_5p3[n5p3] += 1
      d_5p4[n5p4] += 1
      d_5p5[n5p5] += 1

      d_3p1[n3p1] += 1
      d_3p2[n3p2] += 1
      d_3p3[n3p3] += 1
      d_3p4[n3p4] += 1
      d_3p5[n3p5] += 1

      d_md1[mid1] += 1
      d_md2[mid2] += 1
      d_md3[mid3] += 1
      d_md4[mid4] += 1
      d_md5[mid5] += 1

  total_frq_dicercall = sum(d.values())
  ap, cp, gp, tp, gc_content = cal_acgt_percent (total_frq_dicercall, d['A'], d['C'], d['G'], d['T'])
  plenty_features = [total_frq_dicercall, ap, cp, gp, tp, gc_content] 

  many = d_5p1, d_5p2, d_5p3, d_5p4, d_5p5, \
	     d_3p1, d_3p2, d_3p3, d_3p4, d_3p5, \
	     d_md1, d_md2, d_md3, d_md4, d_md5
  for d in many:
    d = normalize (d)
    for k, v in sorted(d.items()):
      plenty_features.append(v)

  return [str(x) for x in plenty_features]

def normalize (d):
  '''
  return the percentage of A,C,G,T of a position.
  '''
  total = sum(d.values())
  d2 = {}
  for k, v in sorted(d.items()):
    if total == 0: d2[k] = 0
    else: d2[k] = round(100 * v/total, 2)
  return d2
  
def get_samdata_from_file (infile):
  #infile = '../input/CONTIG_3__5860000_5865250.sam'
  with open (infile, 'r') as fh:
    samdata = [x.rstrip('\n').split('\t') for x in fh.readlines()]
  return samdata

def title_of_66_features ():
  t = '5p1, 5p2, 5p3, 5p4, 5p5, 3p1, 3p2, 3p3, 3p4, 3p5, md1, md2, md3, md4, md5'.split(', ')
  x = 'A, C, G, T'.split(', ')
  title = 'total_frq_DicerCall, A%, C%, G%, T%, GC%'.split(', ')
  for i in t:
    title += [i + j for j in x]
  return title


def print_argonaute (data, summaryFile):
  with open(summaryFile, 'r') as fh:
    first = fh.readline().rstrip('\n')
        
  outfile = new_summaryFile = summaryFile[:-3] + 'argonautestat.tsv'
  fho = open (outfile, 'w')
  title = title_of_66_features()
  print(first + '\t' + '\t'.join(title), file=fho)
  for i in data:
    print('\t'.join(i), file=fho)
  fho.close()
  return outfile

def demo_samtools_subproc ():
  ''' interval=2:16537288-16538277 '''
  import subprocess as sbp
  DicerCall = 21
  bamfile = '/home/cjwu/scratch/pipeline240320/workdir/sorted_GSM1087987_C0FCSb.bam'
  ch, lb, rb = 2, 16537288, 16538277
  interval = str(ch) + ':' + str(lb) + '-' + str(rb)
  cmd = ['samtools', 'view', bamfile, interval]
  stdout = sbp.check_output(cmd)
  samdata = [x.split() for x in stdout.decode().rstrip().split('\n')]
  plenty_features = argonaute_preference (samdata, DicerCall)
  print(interval, plenty_features)




class argonaute_class ():

  def __init__(self, infile, DicerCall, bamfile):
    self.DicerCall = DicerCall
    self.bamfile = bamfile
    
    with open (infile, 'r') as fh:
        L = fh.readline().rstrip('\n').split('\t')
    ich, ilb, irb = L.index('chr'), L.index('L_bound'), L.index('R_bound')
    self.ich = ich
    self.ilb = ilb
    self.irb = irb
 
  def caculate_acgt_content_from_slice (self, e):
    ''' interval=2:16537288-16538277
        requires sorted and index files (.bam and .bam.bai)
    '''
    ch, lb, rb = e[self.ich], e[self.ilb], e[self.irb] #e[29], e[16], e[17]
    interval = str(ch) + ':' + str(lb) + '-' + str(rb)
    cmd = ['samtools', 'view', self.bamfile, interval]
    stdout = sbp.check_output(cmd)
    samdata = [x.split() for x in stdout.decode().rstrip().split('\n')] 
    if samdata == [[]]: return interval, []
    plenty_features = argonaute_preference (samdata, self.DicerCall)
    return e + plenty_features


if __name__ == "__main__":
  #demo_samtools_subproc () #narval only
  pass


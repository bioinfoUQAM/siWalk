'''
Sequence retrieval and RNA secondary structure prediction.

Wraps samtools faidx for genomic sequence extraction and RNAfold (ViennaRNA)
for secondary structure prediction. Crick-strand sequences are reverse-
complemented before being returned.

Author: Chao-Jung Wu
'''
import subprocess as sbp

def getRevComp(seq):
  """Return the reverse complement of *seq* (DNA, uppercase)."""
  intab, outab = "ACGT", "TGCA"
  trantab = str.maketrans(intab, outab)
  seq = seq.upper()
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

def retrieve(CHR, start, end, genome_file, strand):
  """
  Retrieve a genomic sequence via samtools faidx.

  Returns (coordinate_string, sequence). On error returns (coordinate, 'error').
  Crick-strand ('C') sequences are reverse-complemented.
  """
  coordinate = str(CHR) + ':' + str(start) + '-' + str(end)
  cmd = ['samtools', 'faidx', genome_file, coordinate]
  try:
    stdout = sbp.check_output(cmd, stderr=sbp.STDOUT)
    retrieved_seq = ''.join(stdout.decode().split('\n')[1:]).upper()
    if strand == 'C':
      retrieved_seq = getRevComp(retrieved_seq)
    return coordinate, retrieved_seq
  except sbp.CalledProcessError as e:
    print(f"Command failed with exit status {e.returncode}")
    print(f"Error output: {e.output.decode()}")
    return coordinate, 'error'

def run_RNAfold(seq):
  """
  Fold *seq* with RNAfold (ViennaRNA) and return (dot-bracket structure, MFE).

  Accepts DNA or RNA; RNAfold converts T to U internally.
  Falls back to an unstructured string and MFE=1 if output is malformed.
  """
  cmd1 = ('echo', seq)
  cmd2 = ('RNAfold', '--noPS', '--noLP', '--temp= 20')
  ps = sbp.Popen(cmd1, stdout=sbp.PIPE)
  output = sbp.check_output(cmd2, stdin=ps.stdout)
  output = output.decode().rstrip('\n').split('\n')
  if len(output) < 2: return '.' * len(seq), 1
  seq, data = output
  fold, mfe = data.rstrip(')').split(' (')
  return fold, float(mfe)

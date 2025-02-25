import subprocess as sbp

def getRevComp (seq):
  ''' updated: 2024-Apr-26 '''
  intab, outab = "ACGT", "TGCA"
  trantab = str.maketrans(intab, outab)
  seq = seq.upper()
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

def retrieve (CHR, start, end, genome_file, strand):
  coordinate = str(CHR) + ':' + str(start) + '-' + str(end)  
  cmd = ['samtools', 'faidx', genome_file, coordinate]

  try:
    stdout = sbp.check_output(cmd, stderr=sbp.STDOUT)
    retrieved_seq = ''.join(stdout.decode().split('\n')[1:])
    retrieved_seq = retrieved_seq.upper()
    if strand == 'C':
      retrieved_seq = getRevComp (retrieved_seq)
    return coordinate, retrieved_seq

  except sbp.CalledProcessError as e:
    print(f"Command failed with exit status {e.returncode}")
    print(f"Error output: {e.output.decode()}")
    return coordinate, 'error'
    
def run_RNAfold (seq):
  ''' It does not matter seq has U or T, RNAfold translates T to U within its program.
  Example: echo GUGGAGCUCCUAUCAUUCC| RNAfold
  This requires two subprocesses involving pipe '''
  cmd1 = ('echo', seq)
  cmd2 = ('RNAfold','--noPS', '--noLP', '--temp= 20')
  ps = sbp.Popen(cmd1, stdout=sbp.PIPE)
  output = sbp.check_output(cmd2, stdin=ps.stdout) #b'GUGGAGCUCCUAUCAUUCC\n..((....))......... ( -1.53)\n'
  output = output.decode().rstrip('\n').split('\n')
  if len(output) < 2: return '.'*len(seq), 1
  seq, data = output
  fold, mfe = data.rstrip(')').split(' (') #minimun free energy
  return fold, float(mfe)

def test ():
  seq = 'GUGGAGCUCCUAUCAUUCC'
  fold, mfe = run_RNAfold (seq)
  print(fold, mfe) #..((....))......... -1.53

def test2 ():
  seq = 'GUGGAGCUCCUAUCAUUCC'
  seq = getRevComp (seq)
  fold, mfe = run_RNAfold (seq)
  print(fold, mfe) #.........((....)).. -2.0

def test3 ():
  seq = 'AATACAAAGAGTAGACAATCATCAACCTCATGGTTGGATACAACAAGGAGAGTTCTCAAAATAGATTGAACAATAACATAGATAGGATAATACGACGATAAGTCCATACGTTGGAGACAAGAACAAGATAGTACCTTGCCATAAACAAACGAATGGTGAG'
  fold, mfe = run_RNAfold (seq)
  print(fold, mfe)

# test ()
#test2 ()
'''
miRNA trigger search via miRanda.

Searches each segment's precursor sequence against miRBase miRNAs using the
bundled miRanda binary. Identifies single-hit (any-hit) and two-hit miRNA
triggers (Axtell et al. 2006, Cell 127:565-577).

Author: Chao-Jung Wu
Date:   2024-May-07
'''
import subprocess as sbp
import os.path

AlignLen = 18  # minimum alignment length accepted from miRanda output


def run_miranda(srna, srna_name, precursor, precursor_name, output_tmp):
  '''
  Run miRanda for a single sRNA–precursor pair and return the best hit.

  miRanda handles U/T translation internally and searches only the supplied
  strand; feed the reverse complement to search the trans-strand.

  Returns a list of miRanda output fields for the best hit, or 'nohit'.
  '''
  file1_srna = output_tmp + srna_name
  file2_precursor = output_tmp + precursor_name
  if not os.path.exists(file1_srna):
    with open(file1_srna, 'w') as fh:
      print('>' + srna_name + '(Q.len' + str(len(srna)) + ')\n' + srna, file=fh)
  if not os.path.exists(file2_precursor):
    with open(file2_precursor, 'w') as fh:
      print('>' + precursor_name + '(R.len' + str(len(precursor)) + ')\n' + precursor, file=fh)
  cmd = ['../lib/miranda', file1_srna, file2_precursor, '-noenergy']

  try:
    stdout = sbp.check_output(cmd).decode().rstrip().split('\n')
  except sbp.CalledProcessError as e:
    print(f"Error processing files {file1_srna} and {file2_precursor}: {e}")
    return 'nohit'

  result = next((i for i in stdout if i.startswith('>')), None)
  if result is None: return 'nohit'

  results = []
  for i in stdout:
    if i.startswith('>') and not i.startswith('>>'):
      e = i.split()
      if int(e[8]) >= AlignLen:
        results.append(e)
  if results == []: return 'nohit'

  L1 = results[0]
  for L in results[1:]:
    L1 = [L1[i] if L1[i] == L[i] else f"{L1[i]},{L[i]}" for i in range(len(L1))]
  L1[2] = max(L1[2].split(','))
  return L1


def read_fasta(fasta_file):
  """Read a FASTA file into a {header: sequence} dictionary."""
  d = {}
  current_header = None
  current_sequence = ''
  with open(fasta_file, 'r') as fh:
    for line in fh:
      line = line.strip()
      if line.startswith('>'):
        if current_header:
          d[current_header] = current_sequence
        current_header = line[1:]
        current_sequence = ''
      else:
        current_sequence += line
    if current_header:
      d[current_header] = current_sequence
  return d

def getRevComp(seq):
  """Return the reverse complement of *seq* (DNA, uppercase)."""
  intab, outab = "ACGT", "TGCA"
  trantab = str.maketrans(intab, outab)
  seq = seq.upper()
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

def parse_miRanda_results(outfile, d, species, output_tmp, precursor, precursor_name):
  """Search all species miRNAs against *precursor* and return best-hit summary."""
  fho = open(outfile, 'a')
  Best_miR, BestScore, anyhit_data, twohit_data = '', 0, [], []
  for k, v in d.items():
    if k.startswith(species):
      srna_name, srna = k.split()[0], v
      hitTwice = False
      result = run_miranda(srna, srna_name, precursor, precursor_name, output_tmp)
      if result != 'nohit':
        anyhit_data.append(srna_name)
        if len(result[6].split(',')) > 1:
          hitTwice = True
          twohit_data.append(srna_name)
        data = result + [srna, srna_name, str(hitTwice)]
        print('\t'.join(data), file=fho)
        BestScore_curr = float(result[2])
        if BestScore_curr > BestScore:
          Best_miR, BestScore = srna_name, BestScore_curr
  fho.close()
  if len(anyhit_data) == 0: anyhit_data = ['na']
  if len(twohit_data) == 0: twohit_data = ['na']
  return Best_miR, BestScore, ','.join(anyhit_data), ','.join(twohit_data)


def print_trigger(miranda_data, summaryFile):
  """Append miRanda trigger columns to summaryFile and return the new file path."""
  new_summaryFile = summaryFile[:-3] + 'miRNAtrigger.tsv'
  with open(summaryFile, 'r') as fh:
    first = fh.readline().rstrip('\n')
  first += '\t'.join(', Best_miR, BestScore, anyhit, twohit'.split(', '))
  with open(new_summaryFile, 'w') as fh:
    print(first, file=fh)
    for i in miranda_data:
      print('\t'.join([str(x) for x in i]), file=fh)
  return new_summaryFile


class miRanda_class():

  def __init__(self, infile, ref, species, output_tmp, lib_output, inBasename):
    self.outfile    = lib_output + inBasename + '.miranda_scores.tsv'
    self.d          = read_fasta(ref)
    self.species    = species
    self.output_tmp = output_tmp
    with open(infile, 'r') as fh:
      L = fh.readline().rstrip('\n').split('\t')
    self.iname   = L.index('segment')
    self.ipreseq = L.index('precursor')
    with open(self.outfile, 'w') as fho:
      title = 'Seq1,Seq2,BestScore,TotEnergy(kCal/Mol),Q.start,Q.end,R.start,R.end,AlignLen,UNKNOWN%,UNKNOWN%,srna,srna_name,secondhit'.split(',')
      print('\t'.join(title), file=fho)

  def search_trigger(self, e):
    """Search miRNA triggers for segment *e* and append hit columns."""
    precursor, precursor_name = e[self.ipreseq], e[self.iname]
    Best_miR, BestScore, anyhit_data, twohit_data = parse_miRanda_results(
      self.outfile, self.d, self.species, self.output_tmp, precursor, precursor_name)
    e += [Best_miR, BestScore, anyhit_data, twohit_data]
    return e


if __name__ == "__main__":
  pass

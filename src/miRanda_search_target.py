''' 
2024-May-07
Chao-Jung Wu

@two-hit trigger
Axtell, M.J., Jan, C., Rajagopalan, R. and Bartel, D.P., 2006. A two-hit trigger for siRNA biogenesis in plants. Cell, 127(3), pp.565-577.
    Summary: Two hits mean that there are two or more targeting sites of the same miRNA on the siRNA precursor strand.
'''
import subprocess as sbp
import os.path


AlignLen = 18

def run_miranda (srna, srna_name, precursor, precursor_name, output_tmp):
  ''' miRanda takes care of U/T translation.
      miRanda does not search the trans-strand of the precursor, so to study the targeting on the trans-strand, we need to feed the reverse-complement sequence of the precursor.
      @syntax: miranda file1_srna file2_precursor
      @example: ../lib/miranda ../demo/file1.txt ../demo/file2.txt > ../demo/result.txt
  '''
  file1_srna = output_tmp + srna_name
  file2_precursor = output_tmp + precursor_name
  if not os.path.exists(file1_srna):
    with open (file1_srna, 'w') as fh:
      print ('>' + srna_name + '(Q.len' + str(len(srna)) + ')\n' + srna, file=fh)
  if not os.path.exists(file2_precursor):
    with open (file2_precursor, 'w') as fh:
      print ('>' + precursor_name + '(R.len' + str(len(precursor)) + ')\n' + precursor, file=fh)
  cmd = ['../lib/miranda', file1_srna, file2_precursor, '-noenergy'] # not reporting energy to save time
  
  try:
    stdout = sbp.check_output(cmd).decode().rstrip().split('\n')
  except sbp.CalledProcessError as e:
    print(f"Error processing files {file1_srna} and {file2_precursor}: {e}")
    return 'nohit'

  result = next((i for i in stdout if i.startswith('>')), None)
  if result == None: return 'nohit'

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
  maxScore = max(L1[2].split(','))
  L1[2] = maxScore
  return L1





def read_fasta(fasta_file):
    ''' Read fasta and store as a dictionary
        k: name
        v: squence
    '''
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

def getRevComp (seq):
  ''' updated: 2024-Apr-26 '''
  intab, outab = "ACGT", "TGCA"
  trantab = str.maketrans(intab, outab)
  seq = seq.upper()
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

def parse_miRanda_results(outfile, d, species, output_tmp, precursor, precursor_name):
  fho = open (outfile, 'a')
  Best_miR, BestScore, anyhit_data, twohit_data = '', 0, [], []
  for k, v in d.items():
    if k.startswith(species):
      srna_name, srna= k.split()[0], v
      hitTwice= False
      result = run_miranda (srna, srna_name, precursor, precursor_name, output_tmp)
      if result != 'nohit': #['', '', '', '', '', '', '', '', '', '', '']:
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

def direct_processing():
  ''' Estimate 10 hrs runtime in personal computer.
      For testing , process only the first 9 elements. '''
  outfile = '../output/test_miRNAtrigger.tsv'
  ref = '../dbs/mature.fa'
  species = 'ath'
  summaryFile = '../output/GSM1087987_C0FCSb.library_summary.addref.threshold.tsv'
  new_summaryFile = summaryFile[:-3] + 'miRNAtrigger.tsv'
  output_tmp = '../tmp/'

  with open (summaryFile, 'r') as fh:
    summary_data = [x.rstrip('\n').split('\t') for x in fh.readlines()]
  first = summary_data[0] + 'Best_miR, BestScore, anyhit_data, twohit_data'.split(', ')
  fh = open (new_summaryFile, 'w')
  print('\t'.join(first), file=fh)

  title = 'Seq1,Seq2,TotScore,TotEnergy(kCal/Mol),Q.start,Q.end,R.start,R.end,AlignLen,UNKNOWN%,UNKNOWN%,srna,srna_name,secondhit'.split(',')
  with open (outfile, 'w') as fho:
    print('\t'.join(title), file=fho)

  d = read_fasta(ref)
  for e in summary_data[1:10]: # process only the first 9 elements
    precursor, precursor_name = e[35], e[34]
    Best_miR, BestScore, anyhit_data, twohit_data = parse_miRanda_results(outfile, d, species, output_tmp, precursor, precursor_name)  
    e += [Best_miR, BestScore, anyhit_data, twohit_data]
    print('\t'.join([str(x) for x in e]), file=fh)
  fh.close()  

def demo():
  srna_name = 'tae_m182'
  srna = 'ACCUUCAGGAAGGACUGCAUC'
  precursor = 'AGCAAAATGCACACCAGCCATCCATCTGGTTGAAGGAATACATCATCAGTTATAAACTTATAATCGCACACAAGCACCTAAGGGCTCTTGTGGGCCTTGATGGTACCACGTGGCCAATCAGCAATCCAGTGCTGAAAAACTTAAATAGACATGTATAACCTCATGAGAAAAAGGCATGCGGACAAGGATTAAAAGATTAGCGGATCACCTGCTGCCTTCTAGCCACGCCAGTACACGCGAAGGAAACAAATTATAGAAACCTTCAGACAGCAAACGATGGCGGACTTCGCGCTTGGATTGACCAAGACGGTTGTGGAGGGGACGCTGAGCAGGGTCCAGTTGGCAATCGATGAAGACAACAAGCTGAGGGTGACGGCGAAGCAAGACCTGCGGTACATCACAGCCGAGTTCCAGATTATGCAGTCCTTCCTGAAGGTCGTCAACAAGGAGCGTGCCAATAATGAAGTGGTGAAGACCTGGGTGAAGCAGCTCCGAGACCTTGCCTTTGATGTCGAAGACTGTGTTGAGTTTGTTATCCACTTGGACAACAACTCGCCCTTGACCTGGTTGTGGCGCCTGCTGCCGTCCTGCGTGGCACCACCTCGATCTCGGGACCTGGACAAGGCGGTCGCAGAGTTAAAGGAGCTGAAGGCGAGGGTCGAGGATGTAAGCCAGAGGAACACACGCTATAACCTCATCAATGACTCTGGCTCCCAAGCCAAGACCACCATGCTGACCGACCAGTCCTCCATGGCCACCACCAACCCATCAACATTTCACATGCTAACTGAGGTCTGGGAGGCTGCTGGGAAGCGAGACAGCCTGAGTGACCTACAGAAGTTGATCATGGGTGAAGGCAGTGACCTTCAAACGATCTCTGTATGGGGAAGTACGGCAGCTGATCTTGGGACAATGTCCATTTTCAGCATGATGTATGCTGACCCAGAAATCTGCGGAGCATTCAAAAGACGTGCATGGGTAAAGTTGATGCATCCTTTCAAACCTGATGAGTTCCTCAAGAGCTTGCTTACTCAGTTCTACGTGAGCTCTCACCAAGAAAACATAGATGTGGACCATCTCACGAAGGCCGAGCTCATGCAGCAAGTGAAGGCGCACAAGTATCTCATCATTCTAGAAGAAGTATACAGTGTGGTGGACTGGGATGCCATCCGAGTGTGCCTGCCGGACAATGAGAACGGCAGCCGAATCATCGTGTCGACAAAGCAATTAAGAACTGCACTCTTGTGCACAGGAGGGGCATACCAAGTTTCAGAGCTCAGGAGATTCTCTGATGATCAATCTCTTTGTGCCTTTTCCAAGAAGGTAGGGCGTCAGAGTGGCATGCGAGACCTTATTTGGCAAATAAGGTGTCGTGGTGTGATTTCAGTATGGGGGCTTAGCGAT'
  precursor_name = 'myPrecursor'

  title = 'Seq1,Seq2,BestScore,TotEnergy(kCal/Mol),Q.start,Q.end,R.start,R.end,AlignLen,UNKNOWN%,UNKNOWN%'.split(',')
  result = run_miranda (srna, srna_name, precursor, precursor_name, '../tmp/')
  print(title)
  print(result)
  print()
  # for i in range(len(result)): print(title[i], ':\t', result[i])

def demo_read_fasta():
  # Example usage:
  fasta_file_path = '../dbs/mature.fa'
  species = 'ath'
  d = read_fasta(fasta_file_path)
  for k, v in d.items():
    if k.startswith(species):
      print(k)
      print(v)

def demo_class():
  ''' precursor 3:5862143-5862355 True '''
  precursor = 'GAGGGATAGACAAGGTAGGAGAAAATGACTCGAACGAATTAGAGGTAGAGATAGATATCTATTCTATATTGAGAAGAGATAGAATAGAATCTGTAAAACGAGAAAATAATAAAACGTTTAGAAAGAGATGGGGTCTTACAAGGTCAAGAAAAGGCCTTACAAGGTCAAGAAAAAGAGAATAATGAAATGCATCATCTAGATAAAACACAATAA'
  precursor_name = '3:5862143-5862355'
  ref = '../dbs/mature.fa'
  species = 'ath'
  d = read_fasta(ref)
  
  outfile = 'test_miranda.txt'
  
  fho = open (outfile, 'w')
  title = 'Seq1,Seq2,BestScore,TotEnergy(kCal/Mol),Q.start,Q.end,R.start,R.end,AlignLen,UNKNOWN%,UNKNOWN%,srna,srna_name,secondhit'.split(',')
  print('\t'.join(title), file=fho)
  Best_miR, BestScore, anyhit_data, twohit_data = parse_miRanda_results(outfile, d, species, output_tmp, precursor, precursor_name)
  print('\n==== Best_miR, highest_BestScore:', Best_miR, BestScore)
  print('==== miRNA hitting twice:', twohit_data)
  fho.close()
 
def demo_print():
  summaryFile = '../output/GSM1087987_C0FCSb.library_summary.addref.threshold.tsv'
  outfile = print_trigger ('miranda_data', summaryFile)
  print(outfile)





def print_trigger(miranda_data, summaryFile):
  new_summaryFile = summaryFile[:-3] + 'miRNAtrigger.tsv'
  
  with open(summaryFile, 'r') as fh:
    first = fh.readline().rstrip('\n')
  first += '\t'.join(', Best_miR, BestScore, anyhit, twohit'.split(', '))

  with open(new_summaryFile, 'w') as fh:
    print(first, file=fh)
    for i in miranda_data:
      print('\t'.join([str(x) for x in i]), file=fh)

  return new_summaryFile

    

class miRanda_class ():

  def __init__(self, infile, ref, species, output_tmp, lib_output, inBasename):
    #self.outfile = infile[:-3] + 'miranda_scores.tsv'
    self.outfile = lib_output + inBasename + '.miranda_scores.tsv'
    self.d = read_fasta(ref)
    self.species = species
    self.output_tmp = output_tmp

    with open (infile, 'r') as fh:
        L = fh.readline().rstrip('\n').split('\t')
    iname, ipreseq = L.index('segment'), L.index('precursor')
    self.iname = iname
    self.ipreseq = ipreseq

    fho = open (self.outfile, 'w')
    title = 'Seq1,Seq2,BestScore,TotEnergy(kCal/Mol),Q.start,Q.end,R.start,R.end,AlignLen,UNKNOWN%,UNKNOWN%,srna,srna_name,secondhit'.split(',')
    print('\t'.join(title), file=fho)
    fho.close()

  def search_trigger(self, e):
    ''' 34	segment
        35	precursor '''
    precursor, precursor_name = e[self.ipreseq], e[self.iname]
    Best_miR, BestScore, anyhit_data, twohit_data = parse_miRanda_results(self.outfile, self.d, self.species, self.output_tmp, precursor, precursor_name)
    e.append(Best_miR)
    e.append(BestScore)
    e.append(anyhit_data)
    e.append(twohit_data)
    return e


if __name__ == "__main__":
  #demo_read_fasta()
  ###demo_print()
  #demo_class()

  direct_processing()
  #demo()
  pass




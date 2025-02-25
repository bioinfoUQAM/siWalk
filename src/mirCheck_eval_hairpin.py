'''
Chao-Jung Wu
2024-May-23
'''
import subprocess as sbp


def DEMO_where_siR_is():
  pre = 'CATTAATATAAATAAATGGTTTTTTTACGGGACGGGTTTGGCGGGACGGGTTTGGCAGGACGTTACTTAATAACAATTGTAAACT'
  siR = 'GGGACGGGTTTGGCAGGACGT'
  prefold = '.......................(((((((((((..((((.((((.....)))).)))).))))............)))))))..'
  start = dist_5p = pre.find(siR)
  stop = dist_5p + len(siR)
  
  excised_fold = prefold[start: stop]
  print(siR)
  print(pre[start: stop])
  print(excised_fold)


def DEMO_where_siR_is_part2():
  pre = 'CATTAATATAAATAAATGGTTTTTTTACGGGACGGGTTTGGCGGGACGGGTTTGGCAGGACGTTACTTAATAACAATTGTAAACT'
  siR = 'GGGACGGGTTTGGCAGGACGT'
  start = dist_5p = pre.find(siR)
  extract_pre = pre[start-20: start+50-1]
  new_start = extract_pre.find(siR)
  print(new_start) # new_start = 20


def demo_run_mirCheck (folding, start, stop, option='def'):
  ''' option: def or mey
  cmd: perl eval_mircheck.pl "((((((.((((((....).))))).)))))).........." 46 64 def
  stdout: FAILED_MIR_EXTENSION    0       0
  '''
  cmd = ['perl', 'eval_mircheck.pl', folding, str(start), str(stop), option]
  mircheck_conclu, fback_start, fback_stop = sbp.check_output(cmd).decode().rstrip().split('\t')
  print(mircheck_conclu, fback_start, fback_stop)

def call_mirCheck (folding, start, stop, option='def'):
  cmd = ['perl', 'eval_mircheck.pl', folding, str(start), str(stop), option]
  mircheck_conclu, fback_start, fback_stop = sbp.check_output(cmd).decode().rstrip().split('\t')
  return mircheck_conclu, fback_start, fback_stop

  
  






def print_hairpin_conclu (mircheck_data, summaryFile):
  new_summaryFile = summaryFile[:-3] + 'miRcheck.tsv'

  with open(summaryFile, 'r') as fh:
    first = fh.readline().rstrip('\n')
  first += '\t'.join(', mircheck_conclu, fback_start, fback_stop, mircheck_conclu25, fback_start25, fback_stop25, mircheck_conclu52, fback_start52, fback_stop52'.split(', '))

  with open(new_summaryFile, 'w') as fh:
    print(first, file=fh)
    for i in mircheck_data:
      print('\t'.join([str(x) for x in i]), file=fh)

  return new_summaryFile



class mirCheck_class ():

  def __init__(self, infile, project_path):
    self.option = 'def' # mirCheck options: {def, mey}
    self.plcheck = project_path + 'src/eval_mircheck.pl'
    with open (infile, 'r') as fh:
        L = fh.readline().rstrip('\n').split('\t')
    self.iprefold = L.index('prefold')
    self.idist_5p = L.index('dist_5p')
    self.ieff_seq = L.index('eff_seq')  
    self.iprefold25 = L.index('prefold_200_500')
    self.iprefold52 = L.index('prefold_500_200')

  def run_mirCheck (self, e):
    folding, start, siR = e[self.iprefold], int(e[self.idist_5p]), e[self.ieff_seq]
    stop = start + len(siR)
    cmd = ['perl', self.plcheck, folding, str(start), str(stop), self.option]
    mircheck_conclu, fback_start, fback_stop = sbp.check_output(cmd).decode().rstrip().split('\t')
    e.append(mircheck_conclu)
    e.append(fback_start)
    e.append(fback_stop)

    folding, start = e[self.iprefold25], 200
    stop = start + len(siR)
    cmd = ['perl', self.plcheck, folding, str(start), str(stop), self.option]
    mircheck_conclu, fback_start, fback_stop = sbp.check_output(cmd).decode().rstrip().split('\t')
    e.append(mircheck_conclu)
    e.append(fback_start)
    e.append(fback_stop)

    folding, start = e[self.iprefold52], 500
    stop = start + len(siR)
    cmd = ['perl', self.plcheck, folding, str(start), str(stop), self.option]
    mircheck_conclu, fback_start, fback_stop = sbp.check_output(cmd).decode().rstrip().split('\t')
    e.append(mircheck_conclu)
    e.append(fback_start)
    e.append(fback_stop)
    return e

if __name__ == "__main__":
  #DEMO_where_siR_is()
  DEMO_where_siR_is_part2()
  #demo_run_mirCheck ()
  pass

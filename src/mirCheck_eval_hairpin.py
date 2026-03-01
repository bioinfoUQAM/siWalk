'''
miRcheck hairpin evaluation wrapper.

Calls the Perl script eval_mircheck.pl for each segment to evaluate whether
the precursor sequence forms a valid miRNA-like hairpin. Evaluates three
precursor windows: the segment sequence, a [eff-200, eff+500] window, and
a [eff-500, eff+200] window.

Author: Chao-Jung Wu
Date:   2024-May-23
'''
import subprocess as sbp


def call_mirCheck(folding, start, stop, option='def'):
  """Run miRcheck on a dot-bracket folding and return (conclusion, start, stop)."""
  cmd = ['perl', 'eval_mircheck.pl', folding, str(start), str(stop), option]
  mircheck_conclu, fback_start, fback_stop = sbp.check_output(cmd).decode().rstrip().split('\t')
  return mircheck_conclu, fback_start, fback_stop


def print_hairpin_conclu(mircheck_data, summaryFile):
  """Append miRcheck columns to summaryFile and return the new file path."""
  new_summaryFile = summaryFile[:-3] + 'miRcheck.tsv'
  with open(summaryFile, 'r') as fh:
    first = fh.readline().rstrip('\n')
  first += '\t'.join(', mircheck_conclu, fback_start, fback_stop, mircheck_conclu25, fback_start25, fback_stop25, mircheck_conclu52, fback_start52, fback_stop52'.split(', '))
  with open(new_summaryFile, 'w') as fh:
    print(first, file=fh)
    for i in mircheck_data:
      print('\t'.join([str(x) for x in i]), file=fh)
  return new_summaryFile


class mirCheck_class():

  def __init__(self, infile, project_path):
    self.option  = 'def'  # mirCheck options: {def, mey}
    self.plcheck = project_path + 'src/eval_mircheck.pl'
    with open(infile, 'r') as fh:
      L = fh.readline().rstrip('\n').split('\t')
    self.iprefold   = L.index('prefold')
    self.idist_5p   = L.index('dist_5p')
    self.ieff_seq   = L.index('eff_seq')
    self.iprefold25 = L.index('prefold_200_500')
    self.iprefold52 = L.index('prefold_500_200')

  def run_mirCheck(self, e):
    """Evaluate hairpin for all three precursor windows and append results to *e*."""
    folding, start, siR = e[self.iprefold], int(e[self.idist_5p]), e[self.ieff_seq]
    stop = start + len(siR)
    cmd = ['perl', self.plcheck, folding, str(start), str(stop), self.option]
    mircheck_conclu, fback_start, fback_stop = sbp.check_output(cmd).decode().rstrip().split('\t')
    e += [mircheck_conclu, fback_start, fback_stop]

    # [eff-200, eff+500] window; eff_seq starts at position 200
    folding, start = e[self.iprefold25], 200
    stop = start + len(siR)
    cmd = ['perl', self.plcheck, folding, str(start), str(stop), self.option]
    mircheck_conclu, fback_start, fback_stop = sbp.check_output(cmd).decode().rstrip().split('\t')
    e += [mircheck_conclu, fback_start, fback_stop]

    # [eff-500, eff+200] window; eff_seq starts at position 500
    folding, start = e[self.iprefold52], 500
    stop = start + len(siR)
    cmd = ['perl', self.plcheck, folding, str(start), str(stop), self.option]
    mircheck_conclu, fback_start, fback_stop = sbp.check_output(cmd).decode().rstrip().split('\t')
    e += [mircheck_conclu, fback_start, fback_stop]
    return e


if __name__ == "__main__":
  pass

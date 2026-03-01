'''
Per-contig SAM alignment parser and phasing score calculator.

For each 5250-nt contig, scans every expressed position with a cycle window
(default: 9 cycles × 21 nt) and computes phasing statistics (p, k, maxf, etc.)
along with Howell, Chen, and Guo phase scores. Low-signal positions are
filtered before scoring (k < 3 or maxf < 4; contigs with < 10 expressed
positions are skipped).

Author: Chao-Jung Wu
Date:   2024-Mar-09
'''
import sys
import calculate_Howell_Chen_Guo_scores as cs
import retrieve_seq as ret


def accumulate_count(d, k):
  if k not in d.keys(): d[k] = 1
  else: d[k] += 1
  return d

def get_extra_count(pos, d, param):
  '''
  Sum DicerCall-nt frequencies within ±Dicer_relaxation positions of *pos*.

  @pos   : the central phase position.
  @d     : Watson or Crick pos-freq dictionary.
  @param : [cycle, DicerCall, Dicer_relaxation].
  '''
  cycle, DicerCall, Dicer_relaxation = param
  f = Dicer_relaxation
  extra_count = 0
  for i in range(pos - f, pos + f + 1):
    if i == pos: continue
    if i in d.keys(): extra_count += d[i]
  return extra_count

def parse_sam_data(samdata, param):
  """
  Build per-strand positional frequency dictionaries from raw SAM records.

  Returns [d_21nt_watson, d_21nt_crick, expressed_positions_21nt,
           d_allnt_watson, d_allnt_crick, contig_frq, CHR].
  """
  cycle, DicerCall, Dicer_relaxation = param
  d_21nt_watson_pos_freq  = {}
  d_21nt_crick_pos_freq   = {}
  d_allnt_watson_pos_freq = {}
  d_allnt_crick_pos_freq  = {}
  contig_frq = len(samdata)

  for i in samdata:
    FLAG, CHR, POS, lenSEQ = i[1], i[2], int(i[3]), len(i[9])
    if FLAG == '0':
      STRAND = '+'
      accumulate_count(d_allnt_watson_pos_freq, POS)
    elif FLAG == '16':
      STRAND = '-'
      POS = POS + 2  # shift Crick position to the Watson-strand phase register
      accumulate_count(d_allnt_crick_pos_freq, POS)
    else:
      STRAND = 'error'
      print('error, strand not determined, FLAG =', FLAG)

    if lenSEQ != DicerCall: continue
    if STRAND == '+':
      accumulate_count(d_21nt_watson_pos_freq, POS)
    elif STRAND == '-':
      accumulate_count(d_21nt_crick_pos_freq, POS)

  positions_expressed_by_21nt = sorted(list(set(
    list(d_21nt_watson_pos_freq.keys()) + list(d_21nt_crick_pos_freq.keys()))))
  return [d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, positions_expressed_by_21nt,
          d_allnt_watson_pos_freq, d_allnt_crick_pos_freq, contig_frq, CHR]

def get_p_k_maxf(param, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq):
  '''
  Compute phasing statistics for a cycle window starting at *pos*.

  @param : [cycle, DicerCall, Dicer_relaxation].
  @pos   : window start position.

  Returns: (p, k, maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq).
    p           : total DicerCall-nt frequency on phase positions (with drift).
    k           : number of phase positions with at least one DicerCall-nt sRNA.
    maxf        : frequency of the most abundant phased sRNA (with drift).
    pos_of_maxf : 2-nt adjusted position of maxf.
    eff_strand/pos/frq : effector coordinate allowing drift.
  '''
  cycle, DicerCall, Dicer_relaxation = param
  p, k, tmp = 0, 0, {}
  for i in range(pos, pos + DicerCall * (cycle - 1) + 1, DicerCall):
    w, c, flag = 0, 0, 0
    if i in d_21nt_watson_pos_freq: w = d_21nt_watson_pos_freq[i]; p += w; flag = 1
    if i in d_21nt_crick_pos_freq:  c = d_21nt_crick_pos_freq[i];  p += c; flag = 1
    if flag == 1: k += 1; flag = 0
    tmp[i] = w + c
  pos_of_maxf = max(tmp, key=tmp.get)
  w = get_extra_count(pos_of_maxf, d_21nt_watson_pos_freq, param)
  c = get_extra_count(pos_of_maxf, d_21nt_crick_pos_freq, param)
  maxf = tmp[pos_of_maxf] + w + c
  p += w + c
  eff_strand, eff_pos, eff_frq = get_effector_coordinate(
    pos_of_maxf, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation)
  return p, k, maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq

def get_n_u(param, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq,
            d_allnt_watson_pos_freq, d_allnt_crick_pos_freq):
  '''
  Compute read-count statistics across the full cycle window.

  Returns: (n, N, u, U, ratio, exp, watson_freq, crick_freq).
    n      : positions expressed by ≥1 DicerCall-nt sRNA.
    N      : positions expressed by ≥1 any-size sRNA.
    u      : total DicerCall-nt frequency in window.
    U      : total frequency of all sizes in window.
    ratio  : Watson / total frequency ratio.
    exp    : dominant strand ('Watson', 'Crick', or 'both').
  '''
  cycle, DicerCall, Dicer_relaxation = param
  n, N, u, watson_freq, crick_freq = 0, 0, 0, 0, 0
  for i in range(pos, pos + DicerCall * (cycle - 1) + 1):
    if i in d_21nt_watson_pos_freq: u += d_21nt_watson_pos_freq[i]; n += 1
    if i in d_21nt_crick_pos_freq:  u += d_21nt_crick_pos_freq[i];  n += 1
    if i in d_allnt_watson_pos_freq: watson_freq += d_allnt_watson_pos_freq[i]; N += 1
    if i in d_allnt_crick_pos_freq:  crick_freq  += d_allnt_crick_pos_freq[i];  N += 1
  U = watson_freq + crick_freq
  ratio = round((watson_freq + 0.00001) / (U + 0.00001), 2)
  if ratio > 0.8:   exp = 'Watson'
  elif ratio < 0.2: exp = 'Crick'
  else:             exp = 'both'
  return n, N, u, U, ratio, exp, watson_freq, crick_freq

def parse_positions_expressed_by_21nt(param, many, genome_file):
  """Score all expressed positions and return a dict of {pos: feature_list}."""
  cycle, DicerCall, Dicer_relaxation = param
  [d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, positions_expressed_by_21nt,
   d_allnt_watson_pos_freq, d_allnt_crick_pos_freq, contig_frq, CHR] = many

  d = {}
  for pos in positions_expressed_by_21nt:
    frqw = d_21nt_watson_pos_freq.get(pos, 0)
    frqc = d_21nt_crick_pos_freq.get(pos, 0)
    freq = frqw + frqc

    p, k, maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq = get_p_k_maxf(
      param, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
    if k < 3 or maxf < 4: continue

    n, N, u, U, ratio, exp, watson_freq, crick_freq = get_n_u(
      param, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq,
      d_allnt_watson_pos_freq, d_allnt_crick_pos_freq)

    Howell_score = round(cs.Howell_Xia_2013(p, u, k), 2)
    Howell_2007  = round(cs.Howell_2007(p, k), 2)
    pvalue       = round(cs.Chen_Xia_2013(cycle, k, n), 6)
    pvalue_b     = round(cs.Chen_Xia_2013(cycle, k, N), 6)
    Guo_score    = round(cs.Guo(u, k, p, maxf), 2)
    Guo_score_b  = round(cs.Guo(U, k, p, maxf), 2)

    mid_cycle_pos = pos + round(DicerCall * (cycle - 1) / 2)
    left_bound, right_bound, ext_k, phased_positions = get_boundaries(
      param, k, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
    length = 0 if right_bound == 0 else right_bound - left_bound + 1

    pos_of_maxf_updated, eff_strand_updated, eff_pos_updated, eff_frq_updated = get_updated_eff_pos(
      phased_positions, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation)

    _, seq160 = ret.retrieve(CHR, pos, pos + 160 - 1, genome_file, eff_strand)
    if seq160 == 'error': continue
    fold, mfe = ret.run_RNAfold(seq160)

    mylist = [pos, freq, frqw, frqc, k, n, N, p, u, U,
              maxf, pos_of_maxf_updated, eff_strand_updated, eff_pos_updated, eff_frq_updated,
              mid_cycle_pos, ext_k, left_bound, right_bound, length,
              Howell_score, Howell_2007, Guo_score, Guo_score_b, pvalue, pvalue_b,
              watson_freq, crick_freq, ratio, exp, contig_frq, CHR,
              seq160, fold, mfe]
    d[pos] = mylist
  return d

def get_updated_eff_pos(phased_positions, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation):
  """Return the effector coordinate re-computed from the extended phased position set."""
  tmp = {}
  for i in phased_positions:
    w = d_21nt_watson_pos_freq.get(i, 0)
    c = d_21nt_crick_pos_freq.get(i, 0)
    tmp[i] = w + c
  pos_of_maxf_updated = max(tmp, key=tmp.get)
  eff_strand_updated, eff_pos_updated, eff_frq_updated = get_effector_coordinate(
    pos_of_maxf_updated, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation)
  return pos_of_maxf_updated, eff_strand_updated, eff_pos_updated, eff_frq_updated

def get_phase_pos_in_window(pos, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq):
  """Return the phase positions that have at least one read within a cycle window."""
  collect = []
  for next_pos in range(pos, pos + DicerCall * (cycle - 1) + 1, DicerCall):
    if next_pos in d_21nt_watson_pos_freq or next_pos in d_21nt_crick_pos_freq:
      collect.append(next_pos)
  return collect

def get_boundaries(param, k, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq):
  """Extend the segment boundaries beyond the initial window while phasing signal persists."""
  if k < 3: return 0, 0, k, []

  cycle, DicerCall, Dicer_relaxation = param

  collect = get_phase_pos_in_window(pos, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)

  pos_tmp = pos + DicerCall * cycle
  collect_tmp = get_phase_pos_in_window(pos_tmp, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
  while len(collect_tmp) > 1:
    collect += collect_tmp
    pos_tmp += DicerCall * cycle
    collect_tmp = get_phase_pos_in_window(pos_tmp, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)

  pos_tmp = pos - DicerCall * cycle
  collect_tmp = get_phase_pos_in_window(pos_tmp, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
  while len(collect_tmp) > 1:
    collect += collect_tmp
    pos_tmp -= DicerCall * cycle
    collect_tmp = get_phase_pos_in_window(pos_tmp, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)

  right_bound = max(collect) + DicerCall
  left_bound  = max(min(collect) - DicerCall, 0)
  ext_k = len(collect)
  return left_bound, right_bound, ext_k, collect

def print_stat_table(d_pos_phaseStat, outfile):
  """Write the phasing statistics table to *outfile*; skip if fewer than 10 positions."""
  if len(d_pos_phaseStat) < 10: return

  fho = open(outfile, 'w')
  my_col_names  = 'pos, freq, frqw, frqc, k, n, N, p, u, U'
  my_col_names += ', maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq'
  my_col_names += ', mid_cyc, ext_k, L_bound, R_bound, length'
  my_col_names += ', Howell, Howellb, Guo, Guo_b, pval, pval_b'
  my_col_names += ', Wfreq, Cfreq, ratio, dominant_strand, contig_frq, chromosome'
  my_col_names += ', seq_for_scanning_MEF, fold, mfe'
  print('\t'.join([str(x) for x in my_col_names.split(', ')]), file=fho)
  for pos, mylist in d_pos_phaseStat.items():
    print('\t'.join([str(x) for x in mylist]), file=fho, flush=True)
  fho.close()

def get_samdata_from_file(infile):
  """Read a SAM file and return its records as a list of split lines."""
  with open(infile, 'r') as fh:
    samdata = [x.rstrip('\n').split('\t') for x in fh.readlines()]
  return samdata

def get_samdata_from_stdin():
  """Read SAM records from stdin and return as a list of split lines."""
  import sys
  contents = sys.stdin.read().split('\n')
  samdata = [x.split('\t') for x in contents]
  return samdata

def _internal_eff_coor(d, pos, f):
  """Return the position and frequency of the most abundant read within drift range of *pos*."""
  dres, p, q = {}, 0, 0
  for i in range(pos - f, pos + f + 1):
    if i in d: dres[i] = d[i]
  if dres:
    p = max(dres, key=dres.get)
    q = dres[p]
  return p, q

def get_effector_coordinate(pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation):
  """Return (strand, position, frequency) of the most abundant effector sRNA near *pos*."""
  p,  q  = _internal_eff_coor(d_21nt_watson_pos_freq, pos, Dicer_relaxation)
  p2, q2 = _internal_eff_coor(d_21nt_crick_pos_freq,  pos, Dicer_relaxation)
  if q > q2:   s = 'W'
  elif q == q2: s = 'B'
  else:         s, p, q = 'C', p2 - 2, q2
  return s, p, q

def main(param, infile, outfile):
  samdata = get_samdata_from_file(infile)
  many = parse_sam_data(samdata, param)
  d_pos_phaseStat = parse_positions_expressed_by_21nt(param, many, genome_file)
  print_stat_table(d_pos_phaseStat, outfile)


class parse_alignment_class():

  def __init__(self, output_path, param, genome_file):
    self.output_path = output_path
    self.param       = param
    self.genome_file = genome_file

  def caculate_phasing_scores(self, e, inBasename):
    """Compute phasing scores for contig *e* and write a phasingstat file."""
    key, samdata = e[0], e[1]
    many = parse_sam_data(samdata, self.param)
    d_pos_phaseStat = parse_positions_expressed_by_21nt(self.param, many, self.genome_file)
    print_stat_table(d_pos_phaseStat, self.output_path + inBasename + '.' + key + '.phasingstat.tsv')
    return e


if __name__ == '__main__':
  pass

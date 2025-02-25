'''
Chao-Jung Wu
Date: 2024-Mar-09

This script parse the alignment file that was extracted between two positions, hench a contig.
On this contig, examin each expressed position with its 9-cycle window.
Return the statistics of expressed positions after removing a few non-important postions.
Definition of non-importance is to be determined; for example, remove k < 3.


interval="2:16537288-16538277"
samtools view sorted_${filename}.bam $interval > TAS_chr2_16537288_16538277.sam

Returns: p, k, maxf, ...

Verification update: 241114
  - Low potential segments (e.g., k < 3, maxf < 4, u < 10) were filtered before computing the scores.
    -- around line 187: if k < 3 or maxf < 4: continue
    -- around line 273: if len(d_pos_phaseStat) < 10: return # This is u < 10
'''
import sys
import calculate_Howell_Chen_Guo_scores as cs
import retrieve_seq as ret

# processing_huge_genome = True #False

def accumulate_count (d, k):
    if k not in d.keys(): d[k] = 1
    else: d[k] += 1
    return d

def get_extra_count (pos, d, param):
  '''
  f = 0, 1, or 2
  f =2: Get the accumulate 21-nt freqeuncy at positions of k-2, k-1, k+1, k+2
  '''
  cycle, DicerCall, Dicer_relaxation = param
  f = Dicer_relaxation
  extra_count = 0
  for i in range(pos - f, pos + f + 1):
    if i == pos: continue
    if i in d.keys(): extra_count += d[i]
  return extra_count

def parse_sam_data (samdata, param):
  '''
  retrieve information on the 21-nt-only alignments
  '''
  cycle, DicerCall, Dicer_relaxation = param
  d_21nt_watson_pos_freq = {}
  d_21nt_crick_pos_freq = {}
  d_allnt_watson_pos_freq = {}
  d_allnt_crick_pos_freq = {}
  
  contig_frq = len(samdata)
  
  for i in samdata:
    FLAG, CHR, POS, lenSEQ = i[1], i[2], int(i[3]), len(i[9])
    if FLAG == '0': 
      STRAND = '+'
      accumulate_count (d_allnt_watson_pos_freq, POS)
    elif FLAG == '16': 
      STRAND = '-'
      POS = POS + 2 # record the Crick strand freq at the supposed phasing position of Watson strand.
      accumulate_count (d_allnt_crick_pos_freq, POS)
    else: 
      STRAND = 'error'
      print('error, strand not determined, FLAG =', FLAG) 
    
    if lenSEQ != DicerCall: continue
    if STRAND == '+':
      accumulate_count (d_21nt_watson_pos_freq, POS)
    elif STRAND == '-':
      accumulate_count (d_21nt_crick_pos_freq, POS)
  positions_expressed_by_21nt = sorted(list(set(list(d_21nt_watson_pos_freq.keys()) + list(d_21nt_crick_pos_freq.keys()))))
  return [d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, positions_expressed_by_21nt, 
          d_allnt_watson_pos_freq, d_allnt_crick_pos_freq, contig_frq, CHR]

def get_p_k_maxf (param, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq):
  '''
  Calculates various frequencies within a specified range of positions with step DicerCall.

  Parameters:
  - pos: Starting position to begin the calculation.
  - d_21nt_watson_pos_freq: Dictionary containing frequencies for Watson strand positions.
  - d_21nt_crick_pos_freq: Dictionary containing frequencies for Crick strand positions.

  Returns:
  - p: Sum of phased 21-nt frequencies within the specified range.
  - k: Number of keys (positions) within the specified range.
  - maxf: Maximum phased 21-nt frequency value within the specified range, allowing relaxation of the Dicer specificity
  - pos_of_maxf: 2-nt adjusted position
  - eff_pos: real coordinate in the genome, not adjusted for 2-nt

  Note:
  - DicerCall and cycle are assumed to be predefined constants or variables available
    in the context where this function is called.
  '''
  cycle, DicerCall, Dicer_relaxation = param
  p, k, tmp = 0, 0, {}
  # Iterate over the specified range of positions with step DicerCall.
  for i in range(pos, pos + DicerCall*(cycle-1) + 1, DicerCall):
    w, c, flag = 0, 0, 0
    if i in d_21nt_watson_pos_freq.keys():
      w = d_21nt_watson_pos_freq[i]
      p += w
      flag = 1
    if i in d_21nt_crick_pos_freq.keys():
      c = d_21nt_crick_pos_freq[i]
      p += c
      flag = 1
    if flag == 1: k += 1; flag = 0
    tmp[i] = w + c
  pos_of_maxf = max(tmp, key=tmp.get)
  w = get_extra_count (pos_of_maxf, d_21nt_watson_pos_freq, param)
  c = get_extra_count (pos_of_maxf, d_21nt_crick_pos_freq, param)
  maxf = tmp[pos_of_maxf] + w + c # allowing Dicer relaxation
  p += w + c # allowing Dicer relaxation for the position of the most frequency
  
  eff_strand, eff_pos, eff_frq = get_effector_coordinate (pos_of_maxf, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation)
  return p, k, maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq

def get_n_u (param, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, d_allnt_watson_pos_freq, d_allnt_crick_pos_freq):
  '''
  Calculates various counts and ratios based on the frequencies of positions within a specified range.

  Parameters:
  - pos: Starting position to begin the calculation.
  - d_21nt_watson_pos_freq: Dictionary containing frequencies for Watson strand positions (21-nt).
  - d_21nt_crick_pos_freq: Dictionary containing frequencies for Crick strand positions (21-nt).
  - d_allnt_watson_pos_freq: Dictionary containing frequencies for Watson strand positions (all nucleotides).
  - d_allnt_crick_pos_freq: Dictionary containing frequencies for Crick strand positions (all nucleotides).

  Returns:
  - n: Number of keys (positions) expressed by at least one 21-nt sRNA in the 9-cycle window.
  - N: Number of keys (positions) expressed by at least one any-size sRNA in the 9-cycle window.
  - u: Total frequency of 21-nt sRNAs in the 9-cycle window. Also know as total_abundance.
  - ratio: Ratio of frequencies between Watson and Crick strands considering all RNA sizes.
  - exp: Dominated strand expression based on the ratio 
    'one'-strand for a ratio outside 0.2 to 0.8, 'both'-strand otherwise.
    Watson expression is dominated if ratio > 0.8
    Crick expression is dominated if ratio < 0.2

  Note:
  - DicerCall and cycle are assumed to be predefined constants or variables available
    in the context where this function is called.
  '''
  cycle, DicerCall, Dicer_relaxation = param
  n, N, u, watson_freq, crick_freq = 0, 0, 0, 0, 0
  for i in range(pos, pos + DicerCall*(cycle-1) + 1):
    if i in d_21nt_watson_pos_freq.keys():
      e = d_21nt_watson_pos_freq[i]
      u += e
      n += 1
    if i in d_21nt_crick_pos_freq.keys():
      e = d_21nt_crick_pos_freq[i]
      u += e
      n += 1
    if i in d_allnt_watson_pos_freq.keys():
      e = d_allnt_watson_pos_freq[i]
      watson_freq += e
      N += 1
    if i in d_allnt_crick_pos_freq.keys():
      e = d_allnt_crick_pos_freq[i]
      crick_freq += e
      N += 1
  U = watson_freq + crick_freq
  ratio = round( (watson_freq + 0.00001)/(U + 0.00001), 2)
  if ratio > 0.8: exp = 'Watson'
  elif ratio < 0.2: exp = 'Crick'
  else: exp = 'both'
  return n, N, u, U, ratio, exp, watson_freq, crick_freq

def parse_positions_expressed_by_21nt (param, many, genome_file):
  cycle, DicerCall, Dicer_relaxation = param
  [d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, positions_expressed_by_21nt, 
  d_allnt_watson_pos_freq, d_allnt_crick_pos_freq, contig_frq, CHR] = many
    
  d = {}
  for pos in positions_expressed_by_21nt:
    frqw, frqc = 0, 0
    if pos in d_21nt_watson_pos_freq.keys(): frqw = d_21nt_watson_pos_freq[pos]
    if pos in d_21nt_crick_pos_freq.keys(): frqc = d_21nt_crick_pos_freq[pos]
    freq = frqw + frqc
    # if processing_huge_genome:
        # if freq < 10: continue

    p, k, maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq = get_p_k_maxf (param, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
    if k < 3 or maxf < 4: continue

    n, N, u, U, ratio, exp, watson_freq, crick_freq = get_n_u (param, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, d_allnt_watson_pos_freq, d_allnt_crick_pos_freq)
    # if processing_huge_genome:
        # if u < 10: continue

    Howell_score = round(cs.Howell_Xia_2013 (p, u, k), 2)
    Howell_2007 = round(cs.Howell_2007 (p, k), 2)
    pvalue = round(cs.Chen_Xia_2013 (cycle, k, n), 6)
    pvalue_b = round(cs.Chen_Xia_2013 (cycle, k, N), 6)
    Guo_score = round(cs.Guo (u, k, p, maxf), 2)
    Guo_score_b = round(cs.Guo (U, k, p, maxf), 2)

    mid_cycle_pos = pos + round(DicerCall*(cycle-1)/2)
    left_bound, right_bound, ext_k, phased_positions = get_boundaries (param, k, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
    length = 0 if right_bound == 0 else right_bound - left_bound + 1

    pos_of_maxf_updated, eff_strand_updated, eff_pos_updated, eff_frq_updated = get_updated_eff_pos (phased_positions, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation)

    _, seq160 = ret.retrieve (CHR, pos, pos+160-1, genome_file, eff_strand)
    if seq160 == 'error': continue
    fold, mfe = ret.run_RNAfold (seq160)
    
    mylist = [pos, freq, frqw, frqc, k, n, N, p, u, U, 
              #maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq,
              maxf, pos_of_maxf_updated, eff_strand_updated, eff_pos_updated, eff_frq_updated,
              mid_cycle_pos, ext_k, left_bound, right_bound, length, 
              Howell_score, Howell_2007, Guo_score, Guo_score_b, pvalue, pvalue_b, 
              watson_freq, crick_freq, ratio, exp, contig_frq, CHR, 
              seq160, fold, mfe]
    d[pos] = mylist
  return d

def get_updated_eff_pos (phased_positions, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation):
  tmp = {}
  for i in phased_positions:
    w, c, flag = 0, 0, 0
    if i in d_21nt_watson_pos_freq.keys():
      w = d_21nt_watson_pos_freq[i]
      flag = 1
    if i in d_21nt_crick_pos_freq.keys():
      c = d_21nt_crick_pos_freq[i]
      flag = 1
    if flag == 1: flag = 0
    tmp[i] = w + c
  pos_of_maxf_updated = max(tmp, key=tmp.get)
  eff_strand_updated, eff_pos_updated, eff_frq_updated = get_effector_coordinate (pos_of_maxf_updated, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation)
  return pos_of_maxf_updated, eff_strand_updated, eff_pos_updated, eff_frq_updated

def get_phase_pos_in_window (pos, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq):
  collect = []
  for next_pos in range(pos, pos + DicerCall*(cycle-1) + 1, DicerCall):
    if next_pos in d_21nt_watson_pos_freq.keys(): collect.append(next_pos)
    elif next_pos in d_21nt_crick_pos_freq.keys(): collect.append(next_pos)
  return collect

def get_boundaries (param, k, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq):
  ''' Determines boundaries for extracting precursor based on position. 
      Date: 2024-06-19'''
  if k < 3: return 0, 0, k, []

  cycle, DicerCall, Dicer_relaxation = param
  
  collect = get_phase_pos_in_window(pos, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
  pos_tmp = pos + DicerCall*cycle
  collect_tmp = get_phase_pos_in_window(pos_tmp, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
  while len(collect_tmp) > 1: 
    collect += collect_tmp
    pos_tmp += DicerCall*cycle
    collect_tmp = get_phase_pos_in_window(pos_tmp, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)

  pos_tmp = pos - DicerCall*cycle
  collect_tmp = get_phase_pos_in_window(pos_tmp, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
  while len(collect_tmp) > 1: 
    collect += collect_tmp
    pos_tmp -= DicerCall*cycle
    collect_tmp = get_phase_pos_in_window(pos_tmp, DicerCall, cycle, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)

  right_bound = max(collect) + DicerCall
  left_bound = max( min(collect) - DicerCall, 0) #= max() to ensure left_bound > or = 0
  ext_k = len(collect)
  return left_bound, right_bound, ext_k, collect

def print_stat_table (d_pos_phaseStat, outfile):
  if len(d_pos_phaseStat) < 10: return

  fho = open (outfile, 'w')
  my_col_names = 'pos, freq, frqw, frqc, k, n, N, p, u, U'
  my_col_names += ', maxf, pos_of_maxf, eff_strand, eff_pos, eff_frq'
  my_col_names += ', mid_cyc, ext_k, L_bound, R_bound, length'
  my_col_names += ', Howell, Howellb, Guo, Guo_b, pval, pval_b'
  my_col_names += ', Wfreq, Cfreq, ratio, dominant_strand, contig_frq, chromosome'
  my_col_names += ', seq_for_scanning_MEF, fold, mfe'
  print('\t'.join([str(x) for x in my_col_names.split(', ')]), file=fho)
  for pos, mylist, in d_pos_phaseStat.items():
    print('\t'.join([str(x) for x in mylist]), file=fho, flush=True)
  fho.close()

def get_samdata_from_file (infile):
  #infile = '../input/CONTIG_3__5860000_5865250.sam'
  with open (infile, 'r') as fh:
    samdata = [x.rstrip('\n').split('\t') for x in fh.readlines()]
  return samdata
  
def get_samdata_from_stdin ():
  '''
  #Run the following lines in termial:
  infile=../input/CONTIG_3__5860000_5865250.sam; input_data=$(cat $infile)
  printf "%s" "$input_data" | python3 parse_contig_sam.py
  '''
  import sys
  contents = sys.stdin.read().split('\n')
  samdata = [x.split('\t') for x in contents]
  return samdata

def _internal_eff_coor (d, pos, f):
  ''' allow phase drift
  @d: Watson or Crick pos:freq; due to my code design, need to process Watson and Crick separately
  @p: position of the most abundant freq in d
  @q: freq at p in d
  '''
  dres, p, q = {}, 0, 0
  for i in range(pos - f, pos + f + 1):
    if i in d.keys(): dres[i] = d[i]
  if len(dres.keys()) > 0:
    p = max(dres, key=dres.get)
    q = dres[p]
  return p, q

def get_effector_coordinate (pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation):
  ''' true coordinate, suitable to sequence retrival 
  @p: position of the most abundant freq near pos
  @q: freq at p
  @s: strand of the most abundant freq at p
  '''
  p, q = _internal_eff_coor (d_21nt_watson_pos_freq, pos, Dicer_relaxation)
  p2, q2 = _internal_eff_coor (d_21nt_crick_pos_freq, pos, Dicer_relaxation)
  if q > q2: s = 'W'
  elif q == q2: s = 'B'
  else: s, p, q = 'C', p2-2, q2
  return s, p, q

def main(param, infile, outfile):
  samdata = get_samdata_from_file (infile)
  many = parse_sam_data (samdata, param)
  d_pos_phaseStat = parse_positions_expressed_by_21nt (param, many, genome_file)
  print_stat_table (d_pos_phaseStat, outfile)

def UnitTest_main ():
  infile='../UnitTest_parsing_sam_1/CONTIG_3__5860000_5865250.sam'
  outfile='../UnitTest_parsing_sam_1/CONTIG.3__5860000_5865250.phasingstat.tsv'
  cycle = 9; DicerCall = 21; Dicer_relaxation = 2
  param = [cycle, DicerCall, Dicer_relaxation]
  main(param, infile, outfile)

def UnitTest_boundary ():
  infile='../UnitTest_parsing_sam_1/CONTIG_3__5860000_5865250.sam'
  samdata = get_samdata_from_file (infile)

  cycle = 9; DicerCall = 21; Dicer_relaxation = 2
  param = [cycle, DicerCall, Dicer_relaxation]

  many = parse_sam_data (samdata, param)
  [d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, positions_expressed_by_21nt, _,_,_,_] = many

  maxlength = 0
  for pos in positions_expressed_by_21nt:
    _,k,_,_,_,_,_= get_p_k_maxf (param, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
    left_bound, right_bound, ext_k, phased_positions = get_boundaries (param, k, pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq)
    length = 0 if right_bound == 0 else right_bound - left_bound + 1
    print(left_bound, right_bound, length, ext_k, k)
    print()
    if length > maxlength: maxlength = length
  print('line354, maxlength', maxlength)#316

def UnitTest_effector ():
  infile='../UnitTest_parsing_sam_1/CONTIG_3__5860000_5865250.sam'
  samdata = get_samdata_from_file (infile)
  cycle = 9; DicerCall = 21; Dicer_relaxation = 2
  param = [cycle, DicerCall, Dicer_relaxation]
  many = parse_sam_data (samdata, param)
  [d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, _,_,_,_,_] = many
  
  pos = 5862312
  eff_strand, eff_pos, eff_frq = get_effector_coordinate (
                                 pos, d_21nt_watson_pos_freq, d_21nt_crick_pos_freq, Dicer_relaxation)
  print(eff_strand, eff_pos, eff_frq)


class parse_alignment_class ():

  def __init__(self, output_path, param, genome_file):
    self.output_path = output_path
    self.param = param
    self.genome_file = genome_file

  def caculate_phasing_scores (self, e, inBasename):
    key, samdata = e[0], e[1]
    many = parse_sam_data (samdata, self.param)
    d_pos_phaseStat = parse_positions_expressed_by_21nt (self.param, many, self.genome_file)
    print_stat_table (d_pos_phaseStat, self.output_path + inBasename + '.' + key + '.phasingstat.tsv')
    return e








if __name__ == '__main__' :
  #UnitTest_main ()
  UnitTest_boundary ()
  #UnitTest_effector ()
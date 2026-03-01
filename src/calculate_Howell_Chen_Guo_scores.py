'''
Phasing score functions: Howell, Chen, and Guo schemes.

Implements three phasing score methods used to evaluate siRNA production from
genomic loci. All functions operate on statistics computed within a sliding
cycle window (default: 9 cycles × 21 nt).

References:
  Howell et al. (2007) The Plant Cell 19(3):926-42.
  Chen et al. (2007) PNAS 104(9):3318-23.
  Xia et al. (2013) The Plant Cell 25(5):1555-72.
  Guo & Jin (2015) Bioinformatics 31(2):284-6.

Author: Chao-Jung Wu
Date:   2024-Mar-07
'''
import math

def Howell_Xia_2013(p, u, k):
  '''
  Modified Howell scoring scheme from Xia (2013).

  @p : total frequency of DicerCall-nt sRNAs on phase positions in the window,
       allowing drift for the highest peak (±Dicer_relaxation positions).
  @u : total frequency of DicerCall-nt sRNAs on all positions in the window.
  @k : number of phase positions expressed by at least one DicerCall-nt sRNA.
  '''
  if k < 3: return 0
  x = p / (1 + u)
  y = 1 + 10 * x
  z = math.pow(y, k - 2)
  return math.log(z)

def Howell_2007(p, k):
  '''
  Original Howell (2007) scoring scheme.

  @p : total frequency of DicerCall-nt sRNAs on phase positions in the window,
       allowing drift for the highest peak.
  @k : number of phase positions expressed by at least one DicerCall-nt sRNA.
  '''
  if k < 3: return 0
  z = math.pow(1 + p, k - 2)
  return math.log(z)

def Chen_Xia_2013(c, k, n):
  '''
  Modified Chen scoring scheme from Xia (2013); returns a p-value.

  @c : number of cycles in the window.
  @k : number of phase positions expressed by at least one DicerCall-nt sRNA.
  @n : number of all positions expressed by at least one DicerCall-nt sRNA.
  '''
  if k < 3: return 1
  m = c * 2 - 1
  data = []
  for i in range(k, m + 1):
    if n - i < 0: continue
    A = math.comb(20 * m, n - i)
    B = math.comb(m, i)
    C = math.comb(21 * m, n)
    Pr = A * B / C
    data.append(Pr)
  return sum(data)

def Guo(total_abundance, P_number, P_abundance, max_phased_abundance_allowing_error):
  '''
  Modified Guo (2015) scoring scheme.

  @total_abundance                    (u)    : total frequency of DicerCall-nt sRNAs in the window.
  @P_number                           (k)    : number of phase positions with at least one DicerCall-nt sRNA.
  @P_abundance                        (p)    : total frequency of DicerCall-nt sRNAs on phase positions.
  @max_phased_abundance_allowing_error (maxf) : frequency of the most abundant phased sRNA, allowing drift.
  '''
  if P_number < 3: return 0
  P_ratio = max_phased_abundance_allowing_error / total_abundance
  return P_ratio * P_number * math.log(P_abundance)

def demo_three_calculation():
  p, u, k = 26, 20, 7
  phase_score = Howell_Xia_2013(p, u, k)
  print('Howell-phase-score =', phase_score)

  c, k, n = 9, 4, 23
  pvalue = Chen_Xia_2013(c, k, n)
  print('Chen-pvalue =', pvalue)

  u, k, p, max_phased_abundance_allowing_error = 37, 4, 14, 13
  Guo_phase_score = Guo(u, k, p, max_phased_abundance_allowing_error)
  print('Guo-phase-score =', Guo_phase_score)


if __name__ == '__main__':
  demo_three_calculation()

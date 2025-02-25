'''
Date: 2024-March-07

Calculate phase scores using Howell, Chen and Guo schemes.

References:
@ Xia R, Meyers BC, Liu Z, Beers EP, Ye S, Liu Z. MicroRNA superfamilies descended from miR390 and their roles in secondary small interfering RNA biogenesis in eudicots. The Plant Cell. 2013 May 1;25(5):1555-72.
@ Guo Q, Qu X, Jin W. PhaseTank: genome-wide computational identification of phasiRNAs and their regulatory cascades. Bioinformatics. 2015 Jan 15;31(2):284-6.
@ Chen HM, Li YH, Wu SH. Bioinformatic prediction and experimental validation of a microRNA-directed tandem trans-acting siRNA cascade in Arabidopsis. Proceedings of the National Academy of Sciences. 2007 Feb 27;104(9):3318-23.
@ Howell MD, Fahlgren N, Chapman EJ, Cumbie JS, Sullivan CM, Givan SA, Kasschau KD, Carrington JC. Genome-wide analysis of the RNA-DEPENDENT RNA POLYMERASE6/DICER-LIKE4 pathway in Arabidopsis reveals dependency on miRNA-and tasiRNA-directed targeting. The Plant Cell. 2007 Mar 1;19(3):926-42.

#Return the value of 9 raised to the power of 3
# 9^3
base = 9
power= 3
math.pow(base, power)

# Return the natural logarithm
@x    : Value to calculate > 0
@base : Optional. The logarithmic base to use. Default is 'e'
math.log(x, base) or math.log(x)
'''
import math

def Howell_Xia_2013 (p, u, k):
  '''
  This is modified Howell's scoring scheme in Xia's paper.
  
  Suppose we are examining a 9-cycle window; and Dicer call of 21-nt.
  But the size of window and the size of Dicer call do not affect the calculation here.
  
  @p  : total frequency of all 21-nt sRNAs on phase positions in the 9-cycle window, 
         allowing phase drift for the highest peak by adding the abundance of 21-nt sRNA located at shifting 1 or 2 positions.
  @u  : total frequency of all 21-nt sRNAs on all positions (phase and non-phase) in the 9-cycle window. Also know as total_abundance.
  @k  : number of phase positions expressed by at least one 21-nt sRNA.
  '''
  if k < 3: return 0
  x = p/(1+u)
  y = 1 + 10*x
  z = math.pow(y, k-2)
  return math.log(z)

def Howell_2007 (p, k):
  '''
  This is modified Howell's scoring scheme.
  
  Suppose we are examining a 9-cycle window; and Dicer call of 21-nt.
  But the size of window and the size of Dicer call do not affect the calculation here.
  
  @p  : total frequency of all 21-nt sRNAs on phase positions in the 9-cycle window, 
         allowing phase drift for the highest peak by adding the abundance of 21-nt sRNA located at shifting 1 or 2 positions.
  @u  : total frequency of all 21-nt sRNAs on all positions (phase and non-phase) in the 9-cycle window. Also know as total_abundance.
  @k  : number of phase positions expressed by at least one 21-nt sRNA.
  '''
  if k < 3: return 0
  z = math.pow(1+p, k-2)
  return math.log(z)

def Chen_Xia_2013 (c, k, n):
  '''
  This is modified Chen's scoring scheme in Xia's paper.
   
  Suppose we are examining a 9-cycle window; and Dicer call of 21-nt.
  But the size of window and the size of Dicer call do not affect the calculation here.
  
  @c  : number of cycles
  @m  : number of phase positions availble in the 9-cycle window, minus the first position.
        m = cycle*2-1
        For a 9-cycle window, m = 9*2-1 = 17
        For a 8-cycle window, m = 8*2-1 = 15
  @k  : number of phase positions expressed by at least one 21-nt sRNA.
  @n  : number of all positions expressed by at least one 21-nt sRNA in the 9-cycle window.
  '''
  if k < 3: return 1
  m = c*2-1
  data = []
  for i in range(k, m+1): # i = k, k+1, ..., m
    if n - i < 0: continue
    A = math.comb(20*m, n-i)  
    B = math.comb(m, i)  
    C = math.comb(21*m, n)  
    Pr = A * B /C  
    data.append(Pr)
  return sum(data)

def Guo (total_abundance, P_number, P_abundance, max_phased_abundance_allowing_error):
  '''
  This is modified Guo's scoring scheme.
  
  Suppose we are examining a 9-cycle window; and Dicer call of 21-nt.
  But the size of window and the size of Dicer call do not affect the calculation here.
  
  @total_abundance or u  : total frequency of all 21-nt sRNAs on all positions (phase and non-phase) in the 9-cycle window.
  @P_number        or k  : number of phase positions expressed by at least one 21-nt sRNA.
  @P_abundance     or p  : total frequency of all 21-nt sRNAs on phase positions in the 9-cycle window.
  @max_phased_abundance_allowing_error  or maxf: frequency of the most abundant 21-nt sRNA on phase positions in the 9-cycle window, 
            allowing phase drift for the highest peak by adding the abundance of 21-nt sRNA located at shifting 1 or 2 positions.
  '''
  if P_number < 3: return 0
  P_ratio = max_phased_abundance_allowing_error / total_abundance
  return P_ratio * P_number * math.log(P_abundance)

def demo_three_calculation ():
  p, u, k = 26, 20, 7
  phase_score = Howell_Xia_2013 (p, u, k)
  print('Howell-phase-score =', phase_score)

  c, k, n = 9, 4, 23
  #c, k, n = 9, 1, 6
  pvalue = Chen_Xia_2013 (c, k, n)
  print('Chen-pvalue =', pvalue)

  u, k, p, max_phased_abundance_allowing_error = 37, 4, 14, 13
  Guo_phase_score = Guo (u, k, p, max_phased_abundance_allowing_error)
  print('Guo-phase-score =', Guo_phase_score)


if __name__ == '__main__' :
  demo_three_calculation ()
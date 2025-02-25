'''
2024-May-
Chao-Jung Wu
'''
import pandas as pd
import re
from collections import Counter
import os, csv


def mers123 (infile):
  # Check if the input is a string (assume it's a file path)
  if isinstance(infile, str):
    df = pd.read_csv(infile, sep='\t', low_memory=False)
  # If the input is a DataFrame, use it directly
  elif isinstance(infile, pd.DataFrame):
    df = infile
  else:
    raise ValueError("Input must be a file path or a pandas DataFrame")

  df['eff5p1'] = df['eff_seq'].astype(str).str[0]
  df['eff5p2'] = df['eff_seq'].astype(str).str[1]
  df['eff5p3'] = df['eff_seq'].astype(str).str[2]
  df['eff5p4'] = df['eff_seq'].astype(str).str[3]
  df['eff5p5'] = df['eff_seq'].astype(str).str[4]
  df['eff3p1'] = df['eff_seq'].astype(str).str[-1]
  df['eff3p2'] = df['eff_seq'].astype(str).str[-2]
  df['eff3p3'] = df['eff_seq'].astype(str).str[-3]
  df['eff3p4'] = df['eff_seq'].astype(str).str[-4]
  df['eff3p5'] = df['eff_seq'].astype(str).str[-5]
  df['effmd1'] = df['eff_seq'].astype(str).str[9]
  df['effmd2'] = df['eff_seq'].astype(str).str[10]
  df['effmd3'] = df['eff_seq'].astype(str).str[11]
  df['effmd4'] = df['eff_seq'].astype(str).str[12]
  df['effmd5'] = df['eff_seq'].astype(str).str[13]
  #
  df['eff5p1di'] = df['eff_seq'].astype(str).str[0:2]
  df['eff5p2di'] = df['eff_seq'].astype(str).str[1:3]
  df['eff5p3di'] = df['eff_seq'].astype(str).str[2:4]
  df['eff5p4di'] = df['eff_seq'].astype(str).str[3:5]
  df['eff5p5di'] = df['eff_seq'].astype(str).str[4:6]
  df['eff3p1di'] = df['eff_seq'].astype(str).str[-2:]
  df['eff3p2di'] = df['eff_seq'].astype(str).str[-3:-1]
  df['eff3p3di'] = df['eff_seq'].astype(str).str[-4:-2]
  df['eff3p4di'] = df['eff_seq'].astype(str).str[-5:-3]
  df['eff3p5di'] = df['eff_seq'].astype(str).str[-6:-4]
  df['effmd1di'] = df['eff_seq'].astype(str).str[9:11]
  df['effmd2di'] = df['eff_seq'].astype(str).str[10:12]
  df['effmd3di'] = df['eff_seq'].astype(str).str[11:13]
  df['effmd4di'] = df['eff_seq'].astype(str).str[12:14]
  df['effmd5di'] = df['eff_seq'].astype(str).str[13:15]
  #
  df['eff5p1tri'] = df['eff_seq'].astype(str).str[0:3]
  df['eff5p2tri'] = df['eff_seq'].astype(str).str[1:4]
  df['eff5p3tri'] = df['eff_seq'].astype(str).str[2:5]
  df['eff5p4tri'] = df['eff_seq'].astype(str).str[3:6]
  df['eff5p5tri'] = df['eff_seq'].astype(str).str[4:7]
  df['eff3p1tri'] = df['eff_seq'].astype(str).str[-3:]
  df['eff3p2tri'] = df['eff_seq'].astype(str).str[-4:-1]
  df['eff3p3tri'] = df['eff_seq'].astype(str).str[-5:-2]
  df['eff3p4tri'] = df['eff_seq'].astype(str).str[-6:-3]
  df['eff3p5tri'] = df['eff_seq'].astype(str).str[-7:-4]
  df['effmd1tri'] = df['eff_seq'].astype(str).str[9:12]
  df['effmd2tri'] = df['eff_seq'].astype(str).str[10:13]
  df['effmd3tri'] = df['eff_seq'].astype(str).str[11:14]
  df['effmd4tri'] = df['eff_seq'].astype(str).str[12:15]
  df['effmd5tri'] = df['eff_seq'].astype(str).str[13:16]

  df['eff5p1_3p3'] = df['eff_seq'].astype(str).str[0] + df['eff_seq'].astype(str).str[-3]


  if isinstance(infile, str):
    outfile = infile[:-3] + 'mers.tsv'
    df.to_csv(outfile, sep='\t', index=False)
    return outfile
  elif isinstance(infile, pd.DataFrame):
    # print(df.head())
    return df


def get_feature_value_space (df):
  print(f"\nDataset instances NB= {len(df)}")
  print(f"Dataset features NB= {len(df.columns)}")
  print("\nTo save space, only 3 values of each features are shown.")
  print("==== column \tNB= len(unique_values) \tunique_values[:3] \tmaxval \tminval ====")
  for column in df.columns[:10]:
    unique_values = df[column].unique()
    maxval, minval = df[column].max(), df[column].min()
    mylist = [column, len(unique_values), unique_values[:3], maxval, minval]
    print(mylist)

def demo_iter (df):
  df2 = pd.DataFrame()
  for row in df.itertuples():
    print(row)
    for cols in row:
      print(cols)
      break
    break
  for x, y in df.iteritems():
    print(x)
    print(y)
    break

def demo_feature_value_space():
  infile = '../UnitTest_feature_definition/GSM1087987_C0FCSb.library_summary.addref.threshold.argonautestat.miRNAtrigger.tsv'
  outfile = mers123 (infile)

  df = pd.read_csv(outfile, sep='\t')
  get_feature_value_space (df)
  ###demo_iter (df)










def DEMO_where_siR_is():
  pre = 'CATTAATATAAATAAATGGTTTTTTTACGGGACGGGTTTGGCGGGACGGGTTTGGCAGGACGTTACTTAATAACAATTGTAAACT'
  siR = 'GGGACGGGTTTGGCAGGACGT'
  prefold = '.......................(((((((((((..((((.((((.....)))).)))).))))............)))))))..'
  start = dist_5p = pre.find(siR)
  stop = dist_5p + len(siR)
  
  substring_fold = prefold[start: stop]
  print(siR)
  print(pre[start: stop])
  print(substring_fold)
  
  dist_3p = len(pre) - dist_5p - len(siR)



def loop_with_frequency(fold):
    ''' 
    An example loop: (...) 
    @return: number and lengths of loops within active siRNA
    #fold = "(((.((.((((.....)).....))).))"
    '''
    pattern = r'\(\.+\)'
    substrings = re.findall(pattern, fold)
    frequency = Counter(substrings)
    # Initiate fixed size dictionary
    d = {}
    d['All'] = 0
    for i in range(3, 24): d[i] = 0
    # Assign values to existing keys in dictionary
    for pat, frq in frequency.items():
      d['All'] += frq
      k = len(pat)-2
      if k in d.keys(): d[k] += frq
    return d


def bulges_with_frequency(fold):
    ''' 
    An example bulge representation: A(.(   includes all possibilities such as A(.(  A).)  A).(  but not this A(.) because this will be defined in loop.
    @return: number and lengths of bulges within active siRNA

    Note: wrong pattern #pattern = r'[\(\)]\.+[\(\)]' # ['(.(', '(.....)', ').....)', ').)'] #this includes (.....) which is wrong.
    #fold = "(((.((.((((.....)).....))).))"
    '''
    pattern = r'\(\.+\(|\)\.+\)|\)\.+\(' # ['(.(', ').....)', ').)']
    substrings = re.findall(pattern, fold)
    frequency = Counter(substrings)
    # Initiate fixed size dictionary
    d = {}
    d['All'] = 0
    for i in range(1, 24): d[i] = 0
    # Assign values to existing keys in dictionary
    for pat, frq in frequency.items(): 
      d['All'] += frq
      k = len(pat)-2
      if k in d.keys(): d[k] += frq
    return d


   
def lone_pairs(seq, fold):
  '''
  find .(. or .). and associate with the middle nucleotide
  @return: count of T if in moddle of .(.
  '''
  d = {'A':0, 'T':0, 'C':0, 'G':0,}
  patterns = ['.(.', '.).']
  for pattern in patterns:
    matches = find_all_overlapping_matches(fold, pattern)
    for i in matches:
      mid_nucleotide = seq[i+1]
      if mid_nucleotide not in ['A', 'T', 'C', 'G']: continue
      #print(mid_nucleotide + pattern)
      d[mid_nucleotide] += 1
  return d

def lone_bulge(seq, fold):
  '''
  find ).) or (.( or ).( and associate with the middle nucleotide
  @return: count of T if in moddle of .(.
  '''
  d = {'A':0, 'T':0, 'C':0, 'G':0,}
  patterns = [').)', '(.(', ').(']
  for pattern in patterns:
    matches = find_all_overlapping_matches(fold, pattern)
    for i in matches:
      mid_nucleotide = seq[i+1]
      if mid_nucleotide not in ['A', 'T', 'C', 'G']: continue
      #print(mid_nucleotide + pattern)
      d[mid_nucleotide] += 1
  return d
  


def find_all_overlapping_matches(seq, substring):
  ''' 
  # Find all overlapping matches and their starting indices
  seq = 'TCAAGAAAAAGAGAATAATGA'
  substring = "AAA"
  #All occurrences of 'AAA' found at indices: [5, 6, 7] 
  '''
  matches = []
  start = 0
  while start < len(seq):
    start = seq.find(substring, start)
    if start == -1:
        break
    matches.append(start)
    start += 1  # Move to the next character to allow overlapping matches
  #print(f"All occurrences of 'AAA' found at indices: {matches}")
  return matches




# def create_substring_fold():
  # ''' eff_seq, precursor, prefold, dist_5p '''
  # infile = '../UnitTest_feature_definition/GSM1087987_C0FCSb.library_summary.addref.threshold.argonautestat.miRNAtrigger.tsv'
  # df = pd.read_csv(infile, sep='\t')
  # df['stop'] = df['dist_5p'] + df['eff_seq'].apply(len)
  # df['substring_fold'] = df.apply(lambda row: row['prefold'][row['dist_5p']:row['stop']], axis=1)
  # df.drop(['stop'], axis=1, inplace=True)
  # print(df)


def frequency_of_seq_struc_triplet(seq, fold):
  ''' 
  There are 32 combinations. Triplet A, T, C or G times {paired, unpaired} at each of the three positions; so 4 * 2^3 = 32.
  Let's use . to present unpaired, ( to present paired which includes both open and close brackets.
  Examples: AAA(((   AAA((.   AAA(.(  AAA.((   AAA..(   AAA.(.   AAA(..   AAA...
  '''
  # initiate an empty targeted dictionary
  d = {}
  substrings = ['AAA', 'TTT', 'CCC', 'GGG']
  substructs = ['(((', '((.', '(.(', '.((', '..(', '.(.', '(..', '...']
  # Initiate fixed size dictionary
  for substring in substrings:
    for substruct in substructs:
      k = substring + substruct
      d[k] = 0
  
  for substring in substrings:
    matches = find_all_overlapping_matches(seq, substring)
    for match in matches:
      struc = fold[match:match+3].replace(")", "(")
      k = substring + struc
      if k in d.keys(): d[k] += 1
  return d


def number_of_motif(seq, motif):
  ''' motif = AAA '''
  matches = find_all_overlapping_matches(seq, motif)
  return len(matches)




def demo_replace_string():
  mystring = ")(.((..........)))(("
  modified_string = mystring.replace(")", "(")
  print(mystring)
  print(modified_string)


def get_paired_percentage(fold):
  length = len(fold)
  if length == 0: return 0
  dot_count = fold.count('.')
  paired_percentage = (length-dot_count)*100 / length
  return round(paired_percentage, 2)


def get_paired_rolling_average(fold, window):
  if len(fold) == 0: return 0
  count, accumulate = 0, 0
  for i in range(len(fold)-window+1):
    view = fold[i:i+window]
    dot = view.count('.')
    accumulate += (window - dot)
    count += 1
  average = accumulate / count
  return round(average, 2)



def length_largest_bulge(fold):
    # Define the patterns for bulges
    patterns = [r'\(\.+\(', r'\)\.+\)', r'\)\.+\(']
    # Find all matching bulges
    dot_sequences = []
    for pattern in patterns:
        dot_sequences.extend(re.findall(pattern, fold))
    # Extract the dot sequences from the matches
    dot_sequences = [seq[1:-1] for seq in dot_sequences]  # Remove the surrounding brackets
    if dot_sequences:
        # Find the longest bulge by length
        longest_length = max(len(seq) for seq in dot_sequences)
    else:
        longest_length = 0
    return longest_length


def length_longest_loop(fold):
    # This pattern will match sequences of dots enclosed in parentheses
    pattern = re.compile(r'\((\.+)\)')
    # Find all matching loops
    matches = pattern.findall(fold)
    if matches:
        # Find the longest loop by length
        longest_loop = max(matches, key=len)
        longest_length = len(longest_loop)
    else:
        longest_length = 0
    return longest_length

def length_largest_bracket_sequence(fold):
    # Find all sequences of consecutive brackets
    bracket_sequences = re.findall(r'[()]+', fold)
    # Find the longest sequence by length
    if bracket_sequences:
        longest_length = max(len(seq) for seq in bracket_sequences)
    else:
        longest_length = 0
    return longest_length






def main_outtsv (seq, fold):
  ''' Features from Michael Table 1.
  seq =  'TCAAATTTCCCAGGGTAATGAAAAAAAAA'
  fold = "(((.((.((((.....))..)..))).))"  '''

  d1 = frequency_of_seq_struc_triplet(seq, fold)

  d2 = lone_pairs(seq, fold)
  dtmp = {}
  for k, v in d2.items():
    dtmp[k + '.(.'] = v
  d2 = dtmp  

  d3 = lone_bulge(seq, fold)
  dtmp = {}
  for k, v in d3.items():
    dtmp[k + '(.('] = v
  d3 = dtmp
 
  d4 = bulges_with_frequency(fold)
  dtmp = {}
  for k, v in d4.items():
    dtmp['occ_bulge' + str(k) + 'nt'] = v
  d4 = dtmp

  d5 = loop_with_frequency(fold)
  dtmp = {}
  for k, v in d5.items():
    dtmp['occ_loop' + str(k) + 'nt'] = v
  d5 = dtmp

  return (d1, d2, d3, d4, d5)


def check_and_remove_inconsistent_lines(file_path, expected_fields=None, delimiter='\t'):
    """
    import os, csv
    Check the consistency of the number of fields in each line of a TSV file
    and remove lines with inconsistent fields, replacing the original file.

    Parameters:
    - file_path: str, path to the input TSV file
    - expected_fields: int, number of expected fields in each line. If None, inferred from the first line
    - delimiter: str, the delimiter used in the file (default is tab '\t')
    """
    inconsistent_lines = []
    total_lines = 0
    temp_file_path = file_path + 'tmp'
    with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        with open(temp_file_path, 'w', newline='', encoding='utf-8') as temp_file:
            writer = csv.writer(temp_file, delimiter=delimiter)
            for i, row in enumerate(reader):
                total_lines += 1
                if expected_fields is None:
                    expected_fields = len(row)
                if len(row) == expected_fields:
                    writer.writerow(row)
                else:
                    inconsistent_lines.append(i + 1)
                    print('inconsistent_line', total_lines)
    print(f"Filenanme: {file_path}")
    print(f"Expected {expected_fields} fields.")
    print(f"Removed {len(inconsistent_lines)} inconsistent lines out of {total_lines} total lines.")
    os.replace(temp_file_path, file_path)

def MikeTable1_run(infile):
    ''' cite: Michael Leclercq et al, NAR2013 MiRdup. 
    Report pairing and structure statistics at positions -4 to +4 around start and end of the active siRNA (eff_seq)'''
    outfile = infile[:-3] + 'mikeT1.tsv'
    fho = open (outfile, 'w')

    with open (infile, 'r') as fh:
      DATA = [x.rstrip('\n').split('\t') for x in fh.readlines()]
    L = DATA[0]
    ieff_seq, iprec, iprefold = L.index('eff_seq'), L.index('precursor'), L.index('prefold')

    flag, count = 0, 0
    for i in DATA[1:]:
      eff_seq, precursor, prefold = i[ieff_seq], i[iprec], i[iprefold]
      
      start = dist_5p = precursor.find(eff_seq)
      stop = dist_5p + len(eff_seq)

      a = max(start-4, 0)
      b = min(stop+4, len(precursor))
      substring_fold = prefold[a: b]
      ref_seq = precursor[a: b]

      if len(substring_fold) == 0 or start < 0:
        # print(eff_seq)
        # print(precursor)
        # print(start, stop)
        count += 1
      p = get_paired_percentage(substring_fold)

      title = L + ['paired_percentage']
      data = i + [p]

      for window in [3, 5, 7]:
        c = get_paired_rolling_average(substring_fold, window)
        title.append('paired_roll'+ str(window))
        data.append(c)

      h = length_largest_bulge(substring_fold)
      title.append('length_longest_bulge')
      data.append(h)
      
      h = length_longest_loop(substring_fold)
      title.append('length_longest_loop')
      data.append(h)
      
      h = length_largest_bracket_sequence(substring_fold)
      title.append('longest_paired_length')
      data.append(h)
      
      for motif in ['AAA', 'TTT', 'CCC', 'GGG']:
        c = number_of_motif(ref_seq, motif)
        title.append('NBtriplet'+ motif[0])
        data.append(c)

      ds = main_outtsv(ref_seq, substring_fold)
      for d in ds:
        for k, v in d.items():
          if start < 0: data.append(-1)
          else: data.append(v)
          if flag == 0: title.append(k)

      if flag == 0:
        title = '\t'.join([x for x in title])
        print(title, file=fho)
        flag = 1

      data = '\t'.join([str(x) for x in data])
      print(data, file=fho)
      
    fho.close()
    check_and_remove_inconsistent_lines(outfile)
    return outfile

if __name__ == "__main__":
  
  # create_substring_fold()
  #demo_feature_value_space()
  #DEMO_where_siR_is()

  infile = '../UnitTest_feature_definition/GSM1087987_C0FCSb.library_summary.threshold.argonautestat.miRNAtrigger.miRcheck.mers.addref.threshold.tsv'
  outfile = MikeTable1_run(infile)
  print(outfile)

  fold = "(((.((.((((.....))..)..))).))"  
  # get_paired_percentage(fold)

  # window = 3
  # get_paired_rolling_average(fold, window)


  # length_largest_bulge(fold)
  
  # fold = "(((.((.((((.....))..)..))).......))"
  # longest_bulge_length = length_largest_bulge(fold)
  # print(longest_bulge_length)


  # fold = "(((.((.))))))))(((((.....))..)..))).))"
  # longest_bracket_length = length_largest_bracket_sequence(fold)
  # print(longest_bracket_length)
  
  # seq = 'AAAAACTCAAATCTA'
  # motif = 'AAA'
  # number_of_motif(seq, motif)
  pass

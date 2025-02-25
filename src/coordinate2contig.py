def decide_index (p):
    '''
    @p : position
    @return data: [position, index1 (, index2)]
    
    UPDATED: 2024-MAR-01
    Two extremities are both included.
    If two indices are the same, only one is reported.
    data = [p, index1] or [p, index1, index2]
    Note that index1 > index2
    
    index = [0, 5000, 15000, 20000, 25000...]
    end   = [5250, 15250, 20250, 25250, 30250...]
    p = 1:		index = 0
    p = 5000:		index = 0, 5000
    p = 5100:		index = 0, 5000
    p = 5250:		index = 0, 5000
    p = 10001:            index = 5000, 10000
    p = 10250:            index = 5000, 10000
    p = 10251:            index = 10000
    p = 285741846:        index = 285740000
    '''
    n = int(str(p/5000).split('.')[0])
    index1 = 5000*n

    base4 = int(str(p)[-4:])
    if (5000 <= base4 <= 5250) or (0 <= base4 <= 250):
      index2 = index1 - 5000
      if index2 < 0: index2 = 0
    else: index2 = -1

    data = [p, index1]
    if index1 == index2: pass
    elif index2 > -1: data.append(index2)
    return data

def error_1 (data):
  print("ERROR, p = ", data[0], "computed indice are: ", data[1:])

def unit_test ():
  print("Begin unit test")
  
  data = decide_index (1)
  if data == [1, 0]: pass
  else: error_1 (data)

  data = decide_index (5000)
  if data == [5000, 5000, 0]: pass
  else: error_1 (data)
  
  data = decide_index (5100)
  if data == [5100, 5000, 0]: pass
  else: error_1 (data)
  
  data = decide_index (5250)
  if data == [5250, 5000, 0]: pass
  else: error_1 (data)

  data = decide_index (10001)
  if data == [10001, 10000, 5000]: pass
  else: error_1 (data)

  data = decide_index (10250)
  if data == [10250, 10000, 5000]: pass
  else: error_1 (data)

  data = decide_index (10251)
  if data == [10251, 10000]: pass
  else: error_1 (data)

  data = decide_index (285741846)
  if data == [285741846, 285740000]: pass
  else: error_1 (data)
  
  print("End of unit test\n")

def format_contig (ch, index):
  return ch + "__" + str(index) + "_" + str(index + 5250)

def print_contigs (coordinate, data, CONTIG):
  print("COORDINATE: ", coordinate)
  print("[position, index1 (, index2)] = ", data)
  print("This coordinate belogs to:")
  for contig in CONTIG: print("CONTIG: ", contig)
  print()
  
def convert (ch, p, verbose=False):
  ''' @ch : chromosome 
      @p  : position    '''
  coordinate = [ch, p]
  data = decide_index (p)
  CONTIG = []
  for i in data[1:]:
    contig = format_contig (ch, i)
    CONTIG.append(contig)
  if verbose: print_contigs (coordinate, data, CONTIG)
  return CONTIG

def show_some_examples ():
  print("Begin to show some examples")
  print("EXAMPLE ===============")
  convert ("ChrZ", 285741846, verbose=True)
  convert ("ChrZ", 10250, verbose=True)
  print("End of examples\n")

def test ():
  coordinate = [ch, p] = ["chrZ", 5001]
  data = decide_index (p)
  CONTIG = convert (ch, p)
  print_contigs (coordinate, data, CONTIG)

#unit_test ()
#show_some_examples ()
#test ()
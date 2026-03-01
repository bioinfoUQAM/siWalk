'''
Batch effector localization prediction.

Iterates over a TSV of precursor sequences and runs effector localization
prediction (siWalk_predict_siRNA_location) on each entry, writing results
to a timestamped output file. Intended for batch processing without Spark.

Author: Chao-Jung Wu
Date:   2024-Sep-27
'''
import sys
from datetime import datetime
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
import siWalk_predict_siRNA_location as predict

def main_process_batch(DicerCall, infile, tmpdir='../tmp/', species='ath', mirbase_file='../dbs/mature.fa'):
  """Process all precursors in *infile* and write predicted localization to a TSV."""
  outfile = infile[:-3] + timestamp + 'predicted_localization.tsv'
  fho = open(outfile, 'w')
  title = 'segment, eff_seq, precursor, extract121, real_start, real_end, predict_start, predict_end, indication_score'
  print('\t'.join(title.split(', ')), file=fho, flush=True)
  with open(infile, 'r') as fh:
    fh.readline()
    for line in fh:
      name, _, pri, pre, _, _ = data = line.rstrip('\n').split('\t')
      start, end, score, _, _ = predict.run_one_precursor(name, pre, DicerCall, tmpdir, species, mirbase_file)
      data += [start, end, score]
      print('\t'.join([str(x) for x in data]), file=fho, flush=True)
  fho.close()

def run2(tmpdir):
  DicerCall = 21
  infile = '../UnitTest_input1/consistent_segment_seq_precursor__first3instances.tsv'
  main_process_batch(DicerCall, infile, tmpdir)
  print('Processed:', infile)
  print('End of program')


if __name__ == "__main__":
  print('==== usage: python process_batch_localization.py <infile> <tmpdir> ====')
  print('==== example: python process_batch_localization.py dummpy ../tmp/')
  infile, tmpdir = sys.argv[1], sys.argv[2]
  print('infile:', infile)
  print('tmpdir:', tmpdir)

  if infile == "dummpy":
    run2(tmpdir)
  else:
    main_process_batch(21, infile, tmpdir)
  pass

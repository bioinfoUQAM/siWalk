'''
siWalk — Feature generation pipeline (Module A2).

Processes SAM alignment files using PySpark to compute per-contig features
for siRNA precursor prediction. Outputs one *.contig_features.tsv file per
input library.

Author: Chao-Jung Wu
Date:   2024-Mar-20
'''
SOFTWARE_NAME = 'pipeline240320'
import os
import datetime, time

import arg
import coordinate2contig as cc
import parse_alignment_of_a_contig as pa
import summarize_contigs as su
import gc_content_and_argonaute_preference_of_contig as gc
import add_ref_info as addref
import miRanda_search_target as mst
import create_more_features as cmf
import mirCheck_eval_hairpin as hp
import ml_preprocessing as mpp
import ml_onelib as mo



# Helpers for Spark combineByKey aggregation
def to_list(a): return [a]
def append(a, b): a.append(b); return a
def extend(a, b): a.extend(b); return a



def makedirs_reps(reps):
  """Create each directory in *reps* if it does not already exist."""
  for rep in reps:
    if not os.path.exists(rep):
      os.makedirs(rep)

def pyspark_configuration(appName, masterMemory, heartbeat):
  """Initialise and return a SparkContext with tuned memory and timeout settings."""
  from pyspark import SparkConf, SparkContext
  myConf = SparkConf()
  myConf.setAppName(appName)
  myConf.set("spark.driver.maxResultSize", '1500M')  # default 1g
  myConf.set("spark.driver.memory", masterMemory)
  timeout = heartbeat * 12
  myConf.set('spark.network.timeout', str(timeout) + 's')
  myConf.set('spark.executor.heartbeatInterval', str(heartbeat) + 's')
  return SparkContext(conf=myConf)

def copy_file_to_folder(source, folder):
  """Copy *source* file into *folder* and return the destination path."""
  import shutil
  filename = os.path.basename(source)
  destination = folder.rstrip('/') + '/' + filename
  shutil.copy2(source, destination)
  return destination



msg = '==== siWalk ===='

if __name__ == '__main__':

  print('\n' + msg, 'Initiating and verifying parameters ...')
  args, paramDict = arg.init_params()
  print('============================================================')
  for k, v in sorted(paramDict.items()): print(k, ': ', v)
  print('============================================================\n')

  param = [int(args.cycle), int(args.DicerCall), int(args.Dicer_relaxation)]

  partition = int(args.sc_partition)
  sc = pyspark_configuration(SOFTWARE_NAME, args.sc_mstrmemory, int(args.sc_heartbeat))
  sc.setLogLevel('OFF')  # does not fully mute Spark output
  appId = str(args.jobid) or str(sc.applicationId)

  job_output = args.output_path + appId + '/'

  # Register source, database, and library files with the Spark context
  pyfiles  = [f for f in os.listdir(args.project_path + 'src/')
              if os.path.isfile(os.path.join(args.project_path + 'src/', f))
              and (f.endswith('.py') or f.endswith('.pl'))]
  dbsfiles = [f for f in os.listdir(args.project_path + 'dbs/')
              if os.path.isfile(os.path.join(args.project_path + 'dbs/', f))]
  libfiles = [f for f in os.listdir(args.project_path + 'lib/')
              if os.path.isfile(os.path.join(args.project_path + 'lib/', f))]
  for f in pyfiles:  sc.addPyFile(args.project_path + 'src/' + f)
  for f in dbsfiles: sc.addFile(args.project_path + 'dbs/' + f)
  for f in libfiles: sc.addFile(args.project_path + 'lib/' + f)

  makedirs_reps([job_output, args.output_tmp])

  infiles = [f for f in os.listdir(args.input_path)
             if os.path.isfile(os.path.join(args.input_path, f))]
  print('\n' + msg, 'Processing infiles')
  for infile in infiles:
    if infile[-1:] == '~': print(msg + 'omitting infile', infile); continue
    if not infile.endswith('.sam'): print('omitting', infile); continue
    print('\n' + msg, datetime.datetime.now(), '---- Processing alignments (SAM):', infile)

    inBasename = libname = os.path.splitext(infile)[0]
    infile = args.input_path + infile
    lib_output = job_output + inBasename + '/'
    makedirs_reps([lib_output])
    bamfile = args.input_path + 'sorted_' + inBasename + '.bam'
    if not os.path.isfile(bamfile) or not os.path.isfile(bamfile + '.bai'):
      if not os.path.isfile(bamfile + '.csi'):
        print(msg, 'BAM file or index missing, omitting', infile); continue

    align_rdd = sc.textFile("file:///" + infile, partition)
    if align_rdd.isEmpty(): print(infile, 'is empty, skipping'); continue
    # Filter SAM header lines and split tab-delimited fields
    align_rdd = align_rdd.filter(lambda e: not e.startswith('@')).map(lambda e: e.split('\t'))

    # Map each read to its contig key(s); contig format: chr__start_end (e.g. 3__5860000_5865250)
    index_rdd = align_rdd.flatMap(lambda e: ((i, e) for i in cc.convert(e[2], int(e[3]))))

    # Group reads by contig: (contig_key, [read1, read2, ...])
    combineByKey_rdd = index_rdd.combineByKey(to_list, append, extend)

    if float(args.sampling) == -1:
      # Debug mode: restrict to a single test contig
      combineByKey_rdd = combineByKey_rdd.filter(lambda e: e[0] == '3__5860000_5865250')
    elif float(args.sampling) < 1:
      combineByKey_rdd = combineByKey_rdd.sample(False, float(args.sampling), 12345)

    par_obj = pa.parse_alignment_class(lib_output, param, args.genome_file)
    phas_stat_rdd = combineByKey_rdd.map(lambda e: par_obj.caculate_phasing_scores(e, inBasename))

    # Stage 1: phasing statistics
    print(msg, datetime.datetime.now(), '---- Phasing statistics, contigs:', phas_stat_rdd.count())
    summaryFile = su.run_summary(lib_output, args.output_tmp, inBasename, args.genome_file, param[1])

    # Stage 2: nucleotide content and Argonaute preference
    summary_rdd = sc.textFile("file:///" + summaryFile, partition)\
                    .filter(lambda e: not e.startswith('CONTIG')).map(lambda e: e.split('\t'))
    print(msg, datetime.datetime.now(), '---- Argonaute statistics, segments:', summary_rdd.count())
    gc_obj = gc.argonaute_class(summaryFile, param[1], bamfile)
    sliceGC_rdd = summary_rdd.map(lambda e: gc_obj.caculate_acgt_content_from_slice(e))
    summaryFile = gc.print_argonaute(sliceGC_rdd.collect(), summaryFile)

    # Stage 3: miRNA trigger search via miRanda
    gc_summary_rdd = sc.textFile("file:///" + summaryFile, partition)\
                         .filter(lambda e: not e.startswith('CONTIG')).map(lambda e: e.split('\t'))
    print(msg, datetime.datetime.now(), '---- miRNA triggers, segments:', gc_summary_rdd.count())
    mRnd_obj = mst.miRanda_class(summaryFile, args.mirbase_file, args.species, args.output_tmp, lib_output, inBasename)
    miranda_rdd = gc_summary_rdd.map(lambda e: mRnd_obj.search_trigger(e))
    summaryFile = mst.print_trigger(miranda_rdd.collect(), summaryFile)

    # Stage 4: hairpin evaluation via miRcheck
    trgr_summary_rdd = sc.textFile("file:///" + summaryFile, partition)\
                         .filter(lambda e: not e.startswith('CONTIG')).map(lambda e: e.split('\t'))
    print(msg, datetime.datetime.now(), '---- miRcheck hairpin evaluation, segments:', trgr_summary_rdd.count())
    ck_obj = hp.mirCheck_class(summaryFile, args.project_path)
    mircheck_rdd = trgr_summary_rdd.map(lambda e: ck_obj.run_mirCheck(e))
    summaryFile = hp.print_hairpin_conclu(mircheck_rdd.collect(), summaryFile)

    # Stage 5: additional features (no RDD)
    print(msg, datetime.datetime.now(), '---- Adding mer-like features')
    summaryFile = cmf.mers123(summaryFile)

    print(msg, datetime.datetime.now(), '---- Adding MikeTable1 features')
    summaryFile = cmf.MikeTable1_run(summaryFile)

    summaryFile = copy_file_to_folder(summaryFile, lib_output)
    final_summaryFile = lib_output + inBasename + '.contig_features.tsv'
    os.rename(summaryFile, final_summaryFile)

    print('==========================================================')
  print(msg, datetime.datetime.now(), 'End of jobid ' + appId + ' =============')

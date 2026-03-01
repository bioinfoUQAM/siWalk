'''
Argument parsing for the siWalk pipeline.

Auto-detects project_path from the current working directory (must be run from
src/ or workdir/). Provides defaults for all paths and analysis parameters.

Author: Chao-Jung Wu
Date:   2024-Mar-20
'''
import argparse
import os
import sys

SOFTWARE_NAME = 'pipeline240320'

def find_project_path():
  """Return the project root by stripping /src or /workdir from cwd."""
  cwd = os.getcwd()
  project_path = cwd.split('/src')[0].split('/workdir')[0]
  return project_path

def getOpt(parser):
  """Register all CLI arguments and apply default path resolution."""
  project_path = find_project_path()

  parser.add_argument('--sc_partition', default='64')
  parser.add_argument('--sc_mstrmemory', default='20g')
  parser.add_argument('--sc_heartbeat', default='10')
  parser.add_argument('--jobid')
  parser.add_argument('--project_path', default=project_path)
  parser.add_argument('--input_path')
  parser.add_argument('--output_path')
  parser.add_argument('--output_tmp')

  parser.add_argument('--DicerCall', default='21', choices=[str(i) for i in range(18, 31)])
  parser.add_argument('--cycle', default='9', choices=[str(i) for i in range(8, 21)])
  parser.add_argument('--Dicer_relaxation', default='2', choices=[str(i) for i in range(0, 4)])
  parser.add_argument('--genome_file')
  parser.add_argument('--ref_file')
  parser.add_argument('--species', default='ath')
  parser.add_argument('--mirbase_file')

  parser.add_argument('--sampling', default='1', choices=['1', '0.1', '-1'])

  args = parser.parse_args()
  args.project_path = args.project_path.rstrip('/') + '/'
  args.input_path  = args.input_path  or args.project_path + 'input/'
  args.output_path = args.output_path or args.project_path + 'output/'
  args.output_tmp  = args.output_tmp  or args.project_path + 'tmp/'
  args.output_tmp  = args.output_tmp.rstrip('/') + '/'
  args.mirbase_file = args.mirbase_file or args.project_path + 'dbs/mature.fa'

  paramDict = vars(args)
  return args, paramDict


def init_params():
  """Initialise the argument parser and return (args, paramDict)."""
  parser = argparse.ArgumentParser(prog=SOFTWARE_NAME)
  args, paramDict = getOpt(parser)
  return args, paramDict


if __name__ == '__main__':
  args, paramDict = init_params()
  pass

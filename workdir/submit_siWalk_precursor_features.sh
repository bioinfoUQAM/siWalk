#!/bin/bash
#SBATCH --job-name=siWalk_generate_precursor_features
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=110G
#SBATCH --account=yourgroup


## --------------------------------------
## author: Chao-Jung Wu
## --------------------------------------

# Modify yourgroup and the following variables
INPUT_SAM_PATH=/path/to/folder/
genome_file=/path/to/TAIR10_genomic.fa
ref=/path/to/siWalk/dbs/At_phasing_contig_from_literature.info.tsv
species=ath






# Start the program
module load StdEnv/2020
module load viennarna/2.5.1
module load gcc/9.3.0 blast+/2.14.0
module load spark/3.3.0
module load samtools/1.17 bowtie/1.3.0
module load scipy-stack
chmod +x ../lib/miranda
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index statsmodels
pip install --no-index sklearn
export _JAVA_OPTIONS="-Xms3g -Xmx40g"
spark-submit --master local[*] --executor-memory ${SLURM_MEM_PER_NODE}M \
             ../src/pipeline240320.py \
             --jobid $SLURM_JOBID \
             --input_path $INPUT_SAM_PATH \
			 --output_tmp $SLURM_TMPDIR \
			 --genome_file $genome_file \
             --ref_file $ref \
			 --species $species
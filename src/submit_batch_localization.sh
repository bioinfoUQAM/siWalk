#!/bin/bash
#SBATCH --job-name=possiRloc.array
#SBATCH --time=72:00:00 # 00:05:00 => 5 mins # 01:00:00 1hr
#SBATCH --mem=50G
#SBATCH --account=ctb-banire
#SBATCH --mail-user=g39103001@gm.ym.edu.tw
#SBATCH --mail-type=ALL
#SBATCH --output=%j-%x.out
#SBATCH --error=%j-%x.err
#SBATCH --array=1-64

#git checkout submit_batch_localization.sh && git pull origin master && dos2unix submit_batch_localization.sh && sbatch submit_batch_localization.sh


#job-name=possiRloc.mlclassify.HwlRlStrt.array
#input_file=../UnitTest_input1/consistent_segment_seq_precursor.rankedByHowell.filterRealStart_3k.tsv
#prefix=../tmp/consistent_segment_seq_precursor
input_file=../UnitTest_input1/consistent_segment_seq_precursor.rankedByHowell.filterRealStart.tsv
prefix=../tmp/consistent7k_precursor



# # #job-name=negsiRloc.HwlRlStrt.array
# # input_file=../UnitTest_input1/negative_segment_seq_precursor.rankedByHowell.filterRealStart_3k.tsv
# # prefix=../tmp/negative_segment_seq_precursor
# input_file=../UnitTest_input1/negative_segment_seq_precursor.rankedByHowell.filterRealStart.random6386.tsv
# prefix=../tmp/negative7k_precursor


# # #############################################################################
# # # Mask this bloc when doing array
# # # Generate 64 split files
# N=$(wc -l < "$input_file")  # Get total number of lines from input file
# lines_per_file=$((N / 63))  # Calculate lines per file
# header=$(head -n 1 "$input_file")
# tail -n +2 "$input_file" | split -l $lines_per_file -d -a 3 - split_file_
# counter=1
# for file in split_file_*; do
    # new_file="${prefix}_split_${counter}.txt"
    # echo "$header" | cat - "$file" > "$new_file" 
    # rm "$file"
    # counter=$((counter + 1))
# done
# #############################################################################

# Resulting split file as infile for python script
infile=${prefix}_split_$SLURM_ARRAY_TASK_ID.txt

module load StdEnv/2020
module load viennarna/2.5.1 #perl/5.30.2
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
export _JAVA_OPTIONS="-Xms3g -Xmx10g"

tmpdir=$SLURM_TMPDIR/$SLURM_ARRAY_TASK_ID/
mkdir $tmpdir
python process_batch_localization.py $infile $tmpdir

# python process_batch_localization.py dummpy ../tmp/




# #############################################################################
# # Mask this bloc when doing array
# # After processing, merge resulting annotation file
# # Modify
# # key=consistent_segment_seq_precursor
# # key=negative_segment_seq_precursor
# file1=../tmp2/consistent_segment_seq_precursor_split_9.20241013_222348predicted_localization.tsv
# ls $file1


# key=negative7k_precursor
# pattern=../tmp/$key*20241129*predicted_localization.tsv
# ls $pattern | wc -l
# file1=../tmp/negative7k_precursor_split_9.20241129_015749predicted_localization.tsv
# ls $file1
# datetag=20241129_1125



key=consistent7k_precursor
pattern=../tmp/$key*20241129_1*predicted_localization.tsv
ls $pattern | wc -l
file1=../tmp/consistent7k_precursor_split_9.20241129_121107predicted_localization.tsv
ls $file1
datetag=20241129_1129_2






merged=${key}_predicted_localization_$datetag.tsv
ls $merged
head -n 1 $file1 > $merged
tail -n +2 -q $pattern >> $merged &
# #############################################################################











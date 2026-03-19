# siWalk

SiWalk identifies effector siRNA positions within a given precursor sequence using a Random Forest
classifier trained on structural and sequence features. It operates on sequence alone — no expression
data required.

### Citation
Please cite the paper for your work using this tool.

__Title__: Localization of Phased Effector siRNAs within Precursor Sequences Using Feature-Informed Machine Learning

__Cite__: Chao-Jung Wu, Abdoulaye Baniré Diallo

__BibTex__:
```
@article{... ...
}
```



---

## Quick Start: Predict effector siRNA from a precursor sequence

This is the primary functionality of siWalk.
Given a precursor RNA or DNA sequence, siWalk scans it to identify the most likely effector siRNA position.

### 1. Set up the environment

Download siWalk from GitHub: [github/bioinfoUQAM/siWalk](https://github.com/bioinfoUQAM/siWalk).

**Option A — Local machine:**

Install [ViennaRNA](https://www.tbi.univie.ac.at/RNA/#download) (provides `RNAfold`), then:
```
pip install numpy pandas scipy statsmodels scikit-learn matplotlib
chmod +x /path/siWalk/lib/miranda
```

**Option B — HPC (Compute Canada):**
```
cd /path/siWalk/workdir/
module load StdEnv/2020 viennarna/2.5.1 gcc/9.3.0
module load blast+/2.14.0 samtools/1.17 bowtie/1.3.0 scipy-stack
chmod +x ../lib/miranda
virtualenv --no-download env
source env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index statsmodels sklearn
```

### 2. Run the example

The script `siWalk_predict_siRNA_location.py` takes two arguments: a precursor sequence and DicerCall (the expected siRNA length, typically 21).

**Example**
```
cd /path/siWalk/src/
```

```
priseq=CTTGACCTTGTAAGGCCTTTTCTTGACCTTGTAAGACCCCATCTCTTTCTAAACGTTTTATTATTTTCTCGTTTTACAGATTCTATTCTATCTCTTCTCAATATAGAATAGATATCTATCT
DicerCall=21
python siWalk_predict_siRNA_location.py $priseq $DicerCall
```

> **Note**: Both DNA and RNA sequences are accepted. For best performance, limit precursor sequences to ≤ 200 nucleotides.

### 3. Expected output

The predicted position is printed to stdout:
```
yourPrecursor, predicted siRNA start: 75, end: 93, score: -154.733188
```
The stdout shows the top-ranked candidate. The full list of 6 candidates can be found in `yourPrecursor.effector_localization_top6_recommendation.tsv`.

Output files are written to `/path/siWalk/output/<timestamp>/`:

| File | Description |
|------|-------------|
| `yourPrecursor.effector_localization_indication.tsv` | Scores for all candidate positions |
| `yourPrecursor.effector_localization_top6_recommendation.tsv` | Top 6 ranked siRNA candidates |
| `yourPrecursor.effector_localization_indication.png` | Graphical summary of candidate scores |

The top-6 recommendation file has the following columns:

| Column | Description |
|--------|-------------|
| Position p | Candidate start position on the precursor |
| Start S(p) | Start-position score (sum of weighted indications) |
| End E(p) | End-position score |
| Sum of Indications | Combined score (higher is better) |
| Best length | Predicted siRNA length (nt) |
| End position | Predicted end position |
| ML predicted | Whether this position was flagged by the ML classifier |



---

## Dependencies

__Languages__: `python/3.10.2` (scipy, statsmodels, sklearn, numpy)

__Packages__: `spark/3.3.0`, `bowtie/1.3.0`, `samtools/1.17`, `viennarna/2.5.1`, `perl/5.30.1`, `blast+/2.14.0`

__Distributed with siWalk__: `miRanda`, `miRCheck`

__Environments__: The localization mode (Quick Start above) runs on any local machine or login node.
Modules A and X (precursor prediction and model training) require a Spark-enabled server (developed on Compute Canada Narval).



## Supporting Files

| File | Description |
|------|-------------|
| `dbs/demoOnly_background.tsv` | Small subset of the training set for quick testing |
| `dbs/background.tsv` | Full Arabidopsis training set (185,539 segment samples); download from Zenodo (see below) |
| `dbs/TAIR10_genomic.fa` | [Arabidopsis genome](https://www.ebi.ac.uk/ena/browser/view/GCA_000001735.1); download from Zenodo (see below) |
| `dbs/mature.fa` | miRNAs from [miRBase v22.1](https://mirbase.org/download/) |
| `dbs/GSM1087998_C0RC1.txt` | Example sRNA-seq library in readcount format ([GSE44622](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44622)) |
| `dbs/feature_definition.txt` | Definitions of all variables in the analysis |
| `dbs/complete_feature_list.tsv` | Full feature list for reference |
| `model/GBA100.pkl` | Gradient Boosting + AdaBoost model, top 100 all-type features |
| `model/GBAs100.pkl` | Gradient Boosting + AdaBoost model, top 100 structural features |
| `model/RFAs100.pkl` | Random Forest + AdaBoost model, top 100 structural features; download from Zenodo (see below) |
| `model/Arabidopsis_feature_importance_n_correlation.tsv` | Feature evaluation file paired with `GBA100.pkl` |
| `model/Arabidopsis_structure_feature_importance_n_correlation.tsv` | Feature evaluation file paired with `GBAs100.pkl` and `RFAs100.pkl` |

The following large files can be downloaded from Zenodo and placed as indicated:
```bash
curl -L -O "https://zenodo.org/records/17258908/files/TAIR10_genomic.fa?download=1"  # → dbs/TAIR10/
curl -L -O "https://zenodo.org/records/17258908/files/background.tsv?download=1"     # → dbs/
curl -L -O "https://zenodo.org/records/17258908/files/RFAs100.pkl?download=1"        # → model/
```




---

## Auxiliary Modules

The following modules require a Spark-enabled HPC environment.

### Generate the RFAs100 model

Alternatively to downloading, generate it locally (~30 minutes):
```
cd /path/siWalk/src/
annotated_training_file=../dbs/background.tsv
python siWalk_pickle_localization.py $annotated_training_file
```
Rename the resulting file to `RFAs100.pkl` and move it to `/path/siWalk/model/`.



## Module A. Precursor prediction from sRNA-seq libraries

### A1. Create a SAM file

Align an sRNA-seq library to the genome using `bowtie`:
```
genome_file=/path/TAIR10_genomic.fa
index_prefix=/path/bowtie_index/TAIR10_genomic
bowtie-build -f --quiet $genome_file $index_prefix
samtools faidx $genome_file

path_siWalk=/path/siWalk
filename=GSM1087998_C0RC1
infile=$path_siWalk/dbs/$filename.txt   # in readcount format
raw=$path_siWalk/input/$filename.raw
python $path_siWalk/src/rc2raw.py $infile > $raw

samfile=$path_siWalk/input/$filename.sam
bamfile=$path_siWalk/input/sorted_$filename.bam
bowtie -r $raw -x $index_prefix --threads 64 -v 0 --mm -a --sam --no-unal > $samfile
samtools sort $samfile > $bamfile
samtools index $bamfile
```

### A2. Generate features

Edit `submit_siWalk_precursor_features.sh` to set:
```
INPUT_SAM_PATH=/path/siWalk/input
genome_file=/path/TAIR10_genomic.fa
```

Submit the job from the `workdir` folder:
```
cd /path/siWalk/workdir/
sbatch submit_siWalk_precursor_features.sh
```

Output files will be available in `/path/siWalk/output/$jobID/$filename/`.
The main output is the annotation file: `$filename.contig_features.tsv`.

### A3. Predict precursor segments

```
cd /path/siWalk/src/
annotation_file=/path/siWalk/output/$jobID/$filename/$filename.contig_features.tsv
python siWalk_classify_precursors.py $annotation_file
```

The output `prefix.prediction.tsv` includes these key columns:

| Column | Description |
|--------|-------------|
| CONTIG | 5250-nt genomic region (e.g., `3__5860000_5865250`) |
| eff_strand | Strand of the most abundant effector siRNA |
| eff_seq | Sequence of the most abundant effector siRNA |
| segment | Predicted start and end of the siRNA-generating locus |
| Predicted_Class | `True` if predicted as siRNA precursor, `False` otherwise |



## Module X. Train a machine learning model

### X1. Prepare the training dataset

- Compute features for each precursor candidate using **Module A**.
- Add a `consistent` column: `True` for confirmed precursors, `False` for non-precursors.
- The resulting file is the `annotated_training_file` for the next step.

### X2. Generate the model

**Top 100 all-type features** (produces the GBA100-type model, for precursor prediction only):
```
cd /path/siWalk/src/
annotated_training_file=../dbs/background.tsv
python siWalk_pickle_precursor.py $annotated_training_file
```

**Top 100 structural features** (produces the RFAs100-type model, suitable for both precursor and effector prediction):
```
cd /path/siWalk/src/
annotated_training_file=../dbs/background.tsv
python siWalk_pickle_localization.py $annotated_training_file
```

Each run produces:
- `output/prefix_*_model.pkl` — trained model with timestamp
- `output/prefix_*_feature_importance_n_correlation.tsv` — ranked feature importances and correlations

#!/usr/bin/env bash
#SBATCH --job-name=rdp
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4800M
#SBATCH --mail-type=ALL

# set environmental variables to tell various parallel computation libraries
# how many CPU cores we have
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# java is not on the path by default on compute nodes
module load biojava/21

# GNU time is not on the path by default on compute nodes
export PATH="/appl/opt/time/1.9/bin:$PATH"
TIME="$(which time) --verbose"

# define file names
MODEL=rdp_nbc
DATA=../../data
RESULTS=../../results/$MODEL

TRAIN_FILE=$DATA/train_nt_rdp.fasta

# create input files
TAX_FILE=$DATA/train_tax.txt
export MODEL_DIR=train_nt_${SLURM_CPUS_PER_TASK}
mkdir -p $MODEL_DIR

# taxonomy
cut -f2 -d" " $TAX_FILE |
sort |
uniq |
gawk -F'[|]' '
  BEGIN {
    OFS = "*"
    split("kingdom,phylum,class,order,family,genus,species", ranknames, ",")
    nranks = length(ranknames)
  }
  {
    for (i = 1; i <= NF; ++i) {
      if (i == 1) {
        my_parent = "root"
        my_label = $1
      } else {
        my_parent = my_label
        my_label = $i
      }
      parent[NR][i] = my_parent
      label[NR][i] = my_label
    }
    ntaxa = NR
  }
  END {
    id["root"] = 0
    print 0, "root", 0, 0, "root"
    n = 1
    for (r = 1; r <= nranks; r++) {
      for (i = 1; i <= ntaxa; i++) {
        if (! (label[i][r] in id)) {
          id[label[i][r]] = n
          print n, label[i][r], id[parent[i][r]], r, ranknames[r]
          n++
        }
      }
    }
  }' >train.tax

# train the model
$TIME java -Xmx4g\
           -jar dist/classifier.jar\
           train\
           -t train.tax\
           -s $TRAIN_FILE\
           -o $MODEL_DIR

# the properties files needs to be copied from an existing model
unzip -p rdp_classifier_2.14.zip\
      "rdp_classifier_2.14/src/data/classifier/16srrna/rRNAClassifier.properties"\
      >$MODEL_DIR/$MODEL.properties

# classify the test data
for TESTSET in test testshort;
do
  TESTFILE=$DATA/${TESTSET}_nt.fasta
  RAWFILE=${MODEL}_${TESTSET}_${SLURM_CPUS_PER_TASK}.raw
  $TIME java -Xmx4g\
             -jar dist/classifier.jar\
             classify\
             --conf 0\
             -t $MODEL_DIR/$MODEL.properties\
             -o $RAWFILE\
             $TESTFILE
done

# format the test data and write to results directory
mkdir -p $RESULTS
for f in *_${SLURM_CPUS_PER_TASK}.raw;
do
  RESULT_FILE=$RESULTS/${f%.raw}.tsv
  echo "ID	kingdom	Prob_kingdom	phylum	Prob_phylum	class	Prob_class	order	Prob_order	family	Prob_family	genus	Prob_genus	species	Prob_species" >$RESULT_FILE
  cut -f1,6,8,9,11,12,14,15,17,18,20,21,23,24,26 $f >>$RESULT_FILE
done

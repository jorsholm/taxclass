#!/usr/bin/env bash
#SBATCH --job-name=protax-a
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4800M
#SBATCH --cpus-per-task=40
#SBATCH --array=40
#SBATCH --output=protax-a_%a.out
#SBATCH --error=protax-a_%a.out
#SBATCH --mail-type=ALL

# set environmental variables to tell various parallel computation libraries
# how many CPU cores we have
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# GNU time is not on the path by default on compute nodes
export PATH="/appl/opt/time/1.9/bin:$PATH"
TIME="$(which time) --verbose"

# define file names
MODEL=protax-a
DATA=../../data
RESULTS=../../results/$MODEL

# create input files
TAX_FILE=$DATA/train_tax.txt
TRAIN_FILE=$DATA/train_nt_aln_label.fasta
MODEL_DIR=${MODEL}_${SLURM_ARRAY_TASK_ID}

export PROTAX="$(pwd)/scripts"

# taxonomy
cut -f2 -d" " $TAX_FILE |
sort |
uniq |
gawk -F'[|]' '
  {
    for (i = 1; i <= NF; ++i) {
      if (i == 1) {
        my_parent = "root"
        my_label = $1
      } else {
        my_parent = my_label
        my_label = my_label "," $i
      }
      parent[NR][i] = my_parent
      label[NR][i] = my_label
    }
    ntaxa = NR
  }
  END {
    id["root"] = 0
    print 0, 0, 0, "root"
    print 1, 0, 1, "unk"
    n = 2
    lastparent = "root"
    count[0]=1
    for (r = 1; r <= 7; r++) {
      for (i = 1; i <= ntaxa; i++) {
        if (! (label[i][r] in id)) {
          count[r]++
          if (parent[i][r] != lastparent) {
            print n, id[parent[i][r]], r, parent[i][r] ",unk"
            ++n
            lastparent=parent[i][r]
          }
          id[label[i][r]] = n
          print n, id[parent[i][r]], r, label[i][r]
          n++
        }
      }
    }
    for (r = 1; r <= 7; r++) weight[r] = (0.05 * count[7] / count[r-1])
    print weight[1] "," weight[2] "," weight[3] "," weight[4] "," weight[5] "," weight[6] "," weight[7] >"unk_priors"
  }' >$MODEL_DIR/taxonomy

# mapping reference sequences to taxonomy
awk '{print $1, "7", $2}' $TAX_FILE |
tr "|" "," >$MODEL_DIR/seqid2tax

# reference taxonomy
tr "|" "," <$TRAIN_FILE >$MODEL_DIR/refs.aln

# train the model
cd $MODEL_DIR
$TIME ../train_protax.sh

# Some files need to be combined for inference
cat mcmc? >model.param
cat sc? >model.sc
cat rseqs?.numeric >rseqs.numeric

cd ..

# classify the test data
for TESTSET in test testshort;
do
TESTFILE=$DATA/${TESTSET}_nt_aln.fasta
RAWFILE=${MODEL}_${SLURM_ARRAY_TASK_ID}
$time c2/classify_v2 $MODEL_DIR/taxonomy.priors\
                     $MODEL_DIR/refs.aln\
                     $MODEL_DIR/rseqs.numeric\
                     $MODEL_DIR/model.param\
                     $MODEL_DIR/model.sc\
                     0.01\
                     $DATA/test_nt_aln.fasta\
                     >${MODEL}_${TESTSET}_nt_${SLURM_ARRAY_TASK_ID}.raw
done

# format the test data and write to results directory
mkdir -p $RESULTS
for f in *_${SLURM_ARRAY_TASK_ID}.raw;
do
  RESULT_FILE=$RESULTS/${f%.raw}.tsv
  echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	tribe	Prob_tribe	genus	Prob_genus	species	Prob_species" >$RESULT_FILE
  grep -v "seconds" $f |
  sed 's/[^ ]*,//g' |
  tr ' ' '\t' >>$RESULT_FILE
done

#!/usr/bin/env bash
#SBATCH --job-name=protax
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=4800M
#SBATCH --cpus-per-task=1
#SBATCH --array=1
#SBATCH --output=protax_%a.out
#SBATCH --error=protax_%a.out
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

# Add project-specific bin directory, which includes an R installation.
export PATH="/projappl/project_2005718/bin:$PATH"

# define file names
MODEL=protax
DATA=../../data
RESULTS=../../results/$MODEL

RSEQ2TAX_FILE=$DATA/train_protax.tax
TRAIN_FILE=$DATA/train_nt.fasta

# loop for full and trainonly taxonomy
for TAXONOMY in train_tax full_tax;
do
  # create input files
  TAX_FILE=$DATA/$TAXONOMY.txt
  export MODEL_DIR=train_nt_${TAXONOMY}_${SLURM_ARRAY_TASK_ID}
  mkdir -p $MODEL_DIR

  export PROTAX="$(pwd)/protax/protaxscripts"

  # taxonomy
  sed 's/^[^ ]* //' $TAX_FILE |
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
        leafcount[my_parent]++
        if (i == NF) leafcount[my_label]++
      }
      ntaxa = NR
    }
    END {
      id["root"] = 0
      print 0, 0, 0, "root", 1.0
      n = 1
      lastparent = "root"
      for (r = 1; r <= 7; r++) {
        for (i = 1; i <= ntaxa; i++) {
          if (r > length(label[i])) continue
          if (! (label[i][r] in id)) {
            id[label[i][r]] = n
            print n, id[parent[i][r]], r, label[i][r], leafcount[label[i][r]] / ntaxa
            n++
          }
        }
      }
    }' >$MODEL_DIR/taxonomy

  # mapping reference sequences to taxonomy
  cp $RSEQ2TAX_FILE $MODEL_DIR/seqid2tax

  # reference taxonomy
  tr "|" "," <$TRAIN_FILE >$MODEL_DIR/refs.fa

  # train the model
  cd $MODEL_DIR
  $TIME ../train_protax.sh

  cd ..
  
  export MDIR="$(pwd)/$MODEL_DIR"

  # classify the test data
  for TESTSET in test testshort;
  do
    TESTFILE=$DATA/${TESTSET}_nt_aln.fasta
    export ODIR="$MODEL_DIR/${TESTSET}"
    mkdir -p $ODIR
    cp $TESTFILE $ODIR/test.fa
    cd $ODIR
    
    $TIME ../../classify_protax.sh
    cd ../..
  done
done

# format the test data and write to results directory
mkdir -p $RESULTS
for f in *_${SLURM_ARRAY_TASK_ID}.raw;
do
  RESULT_FILE=$RESULTS/${f%.raw}.tsv
  echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	tribe	Prob_tribe	genus	Prob_genus	species	Prob_species" >$RESULT_FILE
  awk -F' ' '
    /seconds$/ {next}
    {
      for (rank = 1; rank <= 7; rank++) {
        best_prob[rank] = 0
        best_taxon[rank] = "NA"
      }
      for (i = 2; i <= NF; i += 2) {
        rank = split($i, a, ",")
        if (rank > 1 && $i !~ best_taxon[rank - 1]) continue
        if ($(i+1) > best_prob[rank]) {
          best_prob[rank] = $(i+1)
          best_taxon[rank] = $i
          for (rank2 = rank + 1; rank2 <= 7; rank2++) {
            best_prob[rank2] = 0
            best_taxon[rank2] = "NA"
          }
        }
      }
      printf "%s", $1
      for (rank = 1; rank <= 7; rank++) {
        sub(/.+,/, "", best_taxon[rank])
        printf "\t%s\t%.6f", best_taxon[rank], best_prob[rank]
      }
      printf "\n"
    }' $f\
   >>$RESULT_FILE
done
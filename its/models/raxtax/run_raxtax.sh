#!/usr/bin/env bash
#SBATCH --job-name=raxtax
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4800M
#SBATCH --mail-type=ALL

# put vsearch on the path
export PATH="$HOME/.cargo/bin:$PATH"

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
MODEL=raxtax
DATA=../../data
RESULTS=../../results/$MODEL

# make sure the results directory exists
mkdir -p $RESULTS

# train the model
TRAIN_ROOT=train_nt_raxtax
TRAIN_FILE=$DATA/$TRAIN_ROOT.fasta
MODEL_DIR=raxtax_$SLURM_CPUS_PER_TASK
MODEL_FILE=$MODEL_DIR/$TRAIN_ROOT.bin

$TIME raxtax --only-db\
             --database-path $TRAIN_FILE\
             --prefix $MODEL_DIR\
             --threads $SLURM_CPUS_PER_TASK\
             --redo

# classify the test sequences
for TEST in test testshort
do
  TEST_ROOT=${TEST}_nt
  TEST_FILE=$DATA/$TEST_ROOT.fasta
  OUTPUT_DIR=$MODEL_DIR/${TEST_ROOT}
  RAW_FILE=$OUTPUT_DIR/raxtax.out
  mkdir -p $OUTPUT_DIR

  $TIME raxtax --database-path $MODEL_FILE\
               --query-file $TEST_FILE\
               --prefix $OUTPUT_DIR\
               --threads $SLURM_CPUS_PER_TASK\
               --redo\
               --pin

  RESULT_FILE=$RESULTS/${MODEL}_${TEST_ROOT}_${SLURM_CPUS_PER_TASK}.tsv

  # write the header
  echo "ID	kingdom	Prob_kingdom	phylum	Prob_phylum	class	Prob_class	"\
       "order	Prob_order	family	Prob_family	genus	Prob_genus	"\
       "species	Prob_species" >$RESULT_FILE
  # write the data
  awk -F'\t' '
    $1 != query {
      if (NR > 1) {
        printf "%s", query
      }
      query = $1
      for (rank = 1; rank <= 7; rank++) {
        if (NR > 1) {
          printf "\t%s\t%.6f", best_taxon[rank], best_prob[rank]
        }
        best_prob[rank] = 0
        best_taxon[rank] = "NA"
      }
      if (NR > 1) {
        printf "\n"
      }
    }
    {
      n_ranks = split($2, taxon, ",")
      n_probs = split($3, prob, ",")
      if (n_ranks != n_probs) {
        print "Warning: number of ranks and probs do not match for query " query " on line " NR >"/dev/stderr"
        exit 1
      }
      for (rank = 1; rank <= n_ranks; rank++) {
        if (prob[rank] > best_prob[rank]) {
          best_taxon[rank] = taxon[rank]
          best_prob[rank] = prob[rank]
        } else if (taxon[rank] !~ best_taxon[rank]) {
          next
        }
      }
    }
  END{
    printf "%s", query
    for (rank = 1; rank <= 7; rank++) {
      printf "\t%s\t%.6f", best_taxon[rank], best_prob[rank]
    }
    printf "\n"
  }' $RAW_FILE >>$RESULT_FILE

done

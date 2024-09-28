#!/usr/bin/env bash
#SBATCH --job-name=sintax
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4800M
#SBATCH --cpus-per-task=40
#SBATCH --array=40
#SBATCH --output=sintax_%a.out
#SBATCH --error=sintax_%a.out
#SBATCH --mail-type=ALL

# put vsearch on the path
export PATH="/projappl/project_2005718/vsearch/bin:$PATH"

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
MODEL=sintax
DATA=../../data
RESULTS=../../results/$MODEL

# make sure the results directory exists
mkdir -p $RESULTS

# train the model
TRAIN_FILE=$DATA/train_nt_sintax.fasta
MODEL_FILE=train_nt_$SLURM_ARRAY_TASK_ID.udb

$TIME vsearch --makeudb_usearch $TRAIN_FILE\
              --output $MODEL_FILE

# classify the test sequences
for TEST in test testshort
do
  TEST_FILE=$DATA/${TEST}_nt.fasta
  RAW_FILE=${MODEL}_${TEST}_nt_$SLURM_ARRAY_TASK_ID.raw

  $TIME vsearch --sintax $TEST_FILE\
                --db $MODEL_FILE\
                --tabbedout $RAW_FILE\
                --threads $SLURM_CPUS_PER_TASK
done

# reformat the output
for f in *.raw
do
  RESULT_FILE=$RESULTS/${f%.raw}.tsv

  # write the header
  echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	"\
       "subfamily	Prob_subfamily	tribe	Prob_tribe	genus	Prob_genus	"\
       "species	Prob_species" >$RESULT_FILE
  # write the data
  sed -r 's/ [^\t]*\t//; s/[kpcofgs]:([^,\t]+)\(([01][.][0-9]+)\),?/\t\1\t\2/g; s/\t\+.*//'\
       $f\
       >>$RESULT_FILE
done

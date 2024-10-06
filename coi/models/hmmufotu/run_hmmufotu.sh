#!/usr/bin/env bash
#SBATCH --job-name=hmmufotu
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --array=40
##SBATCH --gres=gpu:v100:1
#SBATCH --output=hmmufotu_%a.out
#SBATCH --error=hmmufotu_%a.out
#SBATCH --mail-type=ALL

# activate env module for boost
module load boost

# add project executables to PATH
export PATH="/projappl/project_2005718/bin:$PATH"

# gnu time
export PATH="/appl/opt/time/1.9/bin:$PATH"
TIME="$(which time) --verbose"

# define file names
MODEL=hmmufotu_fasttree
DATA=../../data
RESULTS=../../results/$MODEL

# train model (including building the tree
PREFIX=train_nt
TRAIN_FILE=$DATA/${PREFIX}_aln.fasta
TREE_FILE=${PREFIX}.tree
LOG_FILE=${PREFIX}.vft
TAX_FILE=$DATA/train_tax.tsv

$TIME VeryFastTree -nt\
                   -gtr\
                   -gamma\
                   -threads $SLURM_CPUS_PER_TASK\
                   -out $TREE_FILE\
                   -log $LOG_FILE\
                   $TRAIN_FILE

$TIME hmmufotu-build $TRAIN_FILE\
                     $TREE_FILE\
                     -n $PREFIX\
                     -a $TAX_FILE\
                     -s GTR\
                     -V\
                     -k 8\
                     -p $SLURM_CPUS_PER_TASK

# classify test sequences
for test_case in test testshort
do
  TEST_FILE=$DATA/${test_case}_nt.fasta
  RAW_FILE=${MODEL}_$(basename $TEST_FILE .fasta)_$SLURM_ARRAY_TASK_ID.raw

  $TIME mycoai-classify --model $MODEL_FILE\
                      --out $RAW_FILE\
                      --confidence\
                      $TEST_FILE
done

# format results
for f in ${MODEL}_*_$SLURM_ARRAY_TASK_ID.raw
do
  RESULT_FILE=$RESULTS/${f%.raw}.tsv
  echo -n "ID	order	family	subfamily	tribe	genus	" >$RESULT_FILE
  echo -n "species	Prob_order	Prob_family	" >>$RESULT_FILE
  echo "Prob_subfamily	Prob_tribe	Prob_genus	Prob_species" >>$RESULT_FILE
  cut -d","\
      --output-delimiter=$'\t'\
      -f1\
      --complement\
      ${f} |
      tail -n+2\
      >>$RESULT_FILE
done

echo Job $SLURM_JOB_ID complete

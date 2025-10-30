#!/usr/bin/env bash
#SBATCH --job-name=mycoai_bert
#SBATCH --account=project_2005718
#SBATCH --partition=longrun
#SBATCH --time=14-00:00:00
#SBATCH --mem=32G
##SBATCH --gres=gpu:v100:1
#SBATCH --mail-type=ALL

# load pytorch module
module load pytorch/1.12

# activate python virtual environment
source /projappl/project_2010309/mycoai/bin/activate

# gnu time
export PATH="/appl/opt/time/1.9/bin:$PATH"
TIME="$(which time) --verbose"

# define file names
MODEL=mycoai_bert
DATA=../../data
RESULTS=../../results/$MODEL

# train model
TRAIN_FILE=$DATA/train_nt_unite.fasta
MODEL_FILE=train_nt_bert_$SLURM_CPUS_PER_TASK.pt

$TIME mycoai-train --out $MODEL_FILE\
                   $TRAIN_FILE

# classify test sequences
for test_case in test testshort
do
  TEST_FILE=$DATA/${test_case}_nt.fasta
  RAW_FILE=${MODEL}_$(basename $TEST_FILE .fasta)_$SLURM_CPUS_PER_TASK.raw

  $TIME mycoai-classify --model $MODEL_FILE\
                      --out $RAW_FILE\
                      --confidence\
                      $TEST_FILE
done

# format results
for f in ${MODEL}_*_$SLURM_CPUS_PER_TASK.raw
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

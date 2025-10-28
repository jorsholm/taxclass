#!/usr/bin/env bash
#SBATCH --job-name=mycoai_cnn
#SBATCH --account=project_2005718
#SBATCH --partition=gpu
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:v100:1
#SBATCH --output=mycoai_cnn_gpu.out
#SBATCH --error=mycoai_cnn_gpu.out
#SBATCH --mail-type=ALL

# activate python virtual environment
source /projappl/project_2005718/mycoai/bin/activate

# gnu time
export PATH="/appl/opt/time/1.9/bin:$PATH"
TIME="$(which time) --verbose"

# define file names
MODEL=mycoai_cnn
DATA=../../data
RESULTS=../../results/$MODEL
mkdir -p $RESULTS

# train model
TRAIN_FILE=$DATA/train_nt_unite.fasta
MODEL_FILE=train_nt_cnn_gpu.pt

$TIME mycoai-train --out $MODEL_FILE\
                   --base_arch_type CNN\
                   $TRAIN_FILE

# classify full length test sequences
for test_case in test testshort
do
  TEST_FILE=$DATA/${test_case}_nt.fasta
  RAW_FILE=${MODEL}_$(basename $TEST_FILE .fasta)_gpu.raw

  $TIME mycoai-classify --model $MODEL_FILE\
                        --out $RAW_FILE\
                        --confidence\
                        $TEST_FILE
done

# format results
for f in ${MODEL}_*_gpu.raw
do
  RESULT_FILE=$RESULTS/${f%.raw}.tsv
  echo -n "ID	phylum	class	order	family	genus	" >$RESULT_FILE
  echo -n "species	Prob_phylum	Prob_class	Prob_order	" >>$RESULT_FILE
  echo "Prob_family	Prob_genus	Prob_species" >>$RESULT_FILE
  cut  -d","\
       --output-delimiter=$'\t'\
       -f1\
       --complement\
       ${f} |
  tail -n+2\
       >>$RESULT_FILE
done

echo Completed job $SLURM_JOB_ID

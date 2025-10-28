#!/usr/bin/env bash
#SBATCH --job-name=dnabarcoder
#SBATCH --account=project_2005718
#SBATCH --partition=longrun
#SBATCH --time=7-00:00:00
#SBATCH --mem=64G
##SBATCH --gres=gpu:v100:1
#SBATCH --mail-type=ALL

# install:
# python3.9 -m venv /projappl/project_2005718/dnabarcoder
# source /projappl/project_2005718/dnabarcoder/bin/activate
# pip install biopython matplotlib
# git clone https://github.com/vuthuyduong/dnabarcoder.git

# activate python virtual environment
source /projappl/project_2005718/dnabarcoder/bin/activate

# activate blast environment module
module load blast/2.15.0

# set environmental variables to tell various parallel computation libraries
# how many CPU cores we have
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# gnu time
export PATH="/appl/opt/time/1.9/bin:$PATH"
TIME="$(which time) --verbose"

# define file names
MODEL=dnabarcoder
DATA=../../data
RESULTS=../../results/$MODEL

# abbreviate the "executable"
DNABARCODER="python3 -u dnabarcoder/dnabarcoder.py"

# make sure the results directory exists
mkdir -p $RESULTS

#train the nucleotide model
TRAIN_FILE=$DATA/train_nt.fasta
TAX_FILE=$DATA/train_tax.tsv
SIM_DIR=train_nt_$SLURM_CPUS_PER_TASK
SIM_FILE=$SIM_DIR/$(basename ${TRAIN_FILE%.fasta}.sim)

# calculate distances
$TIME $DNABARCODER sim\
  -i $TRAIN_FILE\
  -o $SIM_DIR

superrank=(species genus family order class phylum kingdom)

# predict thresholds
CUTOFF_FILE=${SIM_DIR}/${SIM_DIR}.cutoffs.json
rm -f $CUTOFF_FILE

# first ranks species..phylum
for i in {1..6}
do
  rank=${superrank[0]}
  superrank=( ${superrank[@]:1} )
  $TIME $DNABARCODER predict\
        -i $TRAIN_FILE\
        -c $TAX_FILE\
        -sim $SIM_FILE\
        -st 0.85\
        -et 1.0\
        -s 0.001\
        -prefix $SIM_DIR\
        -o $SIM_DIR\
        -rank $rank

  $TIME $DNABARCODER predict\
        -i $TRAIN_FILE\
        -c $TAX_FILE\
        -sim $SIM_FILE\
        -st 0.85\
        -et 1.0\
        -s 0.001\
        -prefix $SIM_DIR\
        -o $SIM_DIR\
        -rank $rank\
        -higherrank $(IFS=",$IFS"; echo "${superrank[*]}")
done

# kingdom does not have any superranks; also there are only two "kingdoms"
# (actually classes) so we need minGroupNo
$TIME $DNABARCODER predict\
      -i $TRAIN_FILE\
      -c $TAX_FILE\
      -sim $SIM_FILE\
      -st 0.85\
      -et 1.0\
      -s 0.001\
      -minGroupNo 2\
      -prefix $SIM_DIR\
      -o $SIM_DIR\
      -rank "kingdom"

# calculate best thresholds
# this is the final "model"
BEST_CUTOFF_FILE=${CUTOFF_FILE%.json}.best.json
rm -f $BEST_CUTOFF_FILE
$TIME $DNABARCODER best\
      -i $CUTOFF_FILE\
      -c $TAX_FILE\
      -o $SIM_DIR

# classify test sequences
for test_case in test testshort
do
  TEST_FILE=$DATA/${test_case}_nt.fasta
  TEST_PREFIX=${test_case}_nt_$SLURM_CPUS_PER_TASK
  BESTMATCH_FILE=${TEST_PREFIX}.$(basename $TRAIN_FILE .fasta)_BLAST.bestmatch

  $TIME $DNABARCODER search\
        -i $TEST_FILE\
        -r $TRAIN_FILE\
        --prefix $TEST_PREFIX\
        -o .

  $TIME $DNABARCODER classify\
        -i $BESTMATCH_FILE\
        -r $TRAIN_FILE\
        -c $TAX_FILE\
        --cutoffs $BEST_CUTOFF_FILE\
        -o .\
        --prefix $TEST_PREFIX
done

# format the results
for f in *_$SLURM_CPUS_PER_TASK.classified
do
  RESULT_FILE=$RESULTS/${f%.classified}.tsv
  echo "ID	class	order	family	subfamily	tribe	genus	species"\
       >$RESULT_FILE
  cut -d$'\t'\
      -f1,4\
      --output-delimiter=';'\
      $f |
  tail -n+2 |
  sed 's/;[kpcofgs]__/\t/g'\
      >>$RESULT_FILE
done


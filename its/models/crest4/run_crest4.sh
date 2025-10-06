#!/usr/bin/env bash
#SBATCH --job-name=crest4
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4800M
#SBATCH --cpus-per-task=40
#SBATCH --array=40
#SBATCH --output=crest4_%a.out
#SBATCH --error=crest4_%a.out
#SBATCH --mail-type=ALL

set -e

# load crest4 venv
source /projappl/project_2005718/crest4/bin/activate

# load BLAST
module load blast/2.15.0

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
MODEL=crest4
DATA=../../data
RESULTS=../../results/$MODEL

# make sure the results directory exists
mkdir -p $RESULTS

# train the model
TAX_FILE=$DATA/train_tax.txt

declare -A blasttype
blasttype["nt"]="nucl"
blasttype["aa"]="prot"

declare -A blastcmd
blastcmd["nt"]="blastn"
blastcmd["aa"]="blastp"

# classify the test sequences
for TYPE in nt
do
  TRAIN_FILE=$DATA/train_$TYPE.fasta
  MODEL_DIR=train_${TYPE}_$SLURM_ARRAY_TASK_ID

  mkdir -p $MODEL_DIR
  rm -rf $MODEL_DIR/*

  tr '| ' '/\t' <$TAX_FILE |
  sed -r -e 's|\t|\troot/|'\
         -e 's|([^/]*)$|\1\t\1|'\
         >$MODEL_DIR/$MODEL_DIR.tsv
  cp $TRAIN_FILE $MODEL_DIR/$MODEL_DIR.fasta

  cd $MODEL_DIR

  $TIME python3 ../crest4_utils/make_new_crest_db.py $MODEL_DIR.tsv
  $TIME makeblastdb -in $MODEL_DIR.fasta\
                    -out $MODEL_DIR\
                    -dbtype ${blasttype[$TYPE]}

  cd ..


  for TEST in test testshort
  do
    TEST_FILE=$DATA/${TEST}_${TYPE}.fasta
    RAW_DIR=${MODEL}_${TEST}_${TYPE}_$SLURM_ARRAY_TASK_ID
    HITS_FILE=${RAW_DIR}/search.hits
    mkdir -p $RAW_DIR
    rm -rf $RAW_DIR/*
    # run independent blast
    $TIME ${blastcmd[$TYPE]} -query $TEST_FILE\
                             -db $MODEL_DIR/$MODEL_DIR\
                             -num_alignments 100\
                             -outfmt "7 qseqid sseqid bitscore length nident"\
                             -out $HITS_FILE\
                             $([ $TYPE = "nt" ] && echo "-strand plus" || true)\
                             -num_threads $SLURM_CPUS_PER_TASK
    # remove "BLAST processed" lines which mess up BioPython's parser
    sed -ni '/# BLAST processed/!p' $HITS_FILE
    # run crest4 (it always returns 1)
    $TIME crest4 --fasta $TEST_FILE\
                 --search_db $MODEL_DIR/\
                 --search_hits $HITS_FILE\
                 --num_threads $SLURM_CPUS_PER_TASK\
                 --output_dir $RAW_DIR ||
          [ $? = 1 ]
  done
done

# reformat the output
for f in */assignments.txt
do
  RESULT_FILE=$RESULTS/$(dirname $f).tsv

  # write the header
  echo "ID	kingdom	phylum	class	order	family	genus	species" >$RESULT_FILE
  # write the data
  awk -F'[;\t] *'\
      '{
        printf "%s",$1
        for (i=4; i<=10; i++) {
         if (i <= NF) {
          printf "\t%s",$i
         } else {
          printf "\tunclassified"
         }
        }
        print ""
       }'\
       $f\
       >>$RESULT_FILE
done

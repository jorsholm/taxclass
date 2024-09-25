#!/usr/bin/env bash
#SBATCH --job-name=blast
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4800M
#SBATCH --cpus-per-task=4
#SBATCH --array=4
#SBATCH --output=blast_%a.out
#SBATCH --error=blast_%a.out
#SBATCH --mail-type=ALL

# make BLAST available as a module
module load blast

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
MODEL=blast
DATA=../../data
RESULTS=../../results/$MODEL

# make sure the results directory exists
mkdir -p $RESULTS

#train the nucleotide model
TRAIN_FILE=$DATA/train_nt_sintax.fasta
MODEL_FILE=train_nt_$SLURM_ARRAY_TASK_ID

$TIME makeblastdb -in $TRAIN_FILE\
                  -dbtype nucl\
                  -out $MODEL_FILE

# classify the nucleotide test sequences
for TEST in test testshort
do
  TEST_FILE=$DATA/${TEST}_nt.fasta
  RAW_FILE=${MODEL}_${TEST}_nt_$SLURM_ARRAY_TASK_ID.raw

  $TIME blastn -query $TEST_FILE\
               -db $MODEL_FILE\
               -out $RAW_FILE\
               -strand plus\
               -num_threads $SLURM_CPUS_PER_TASK\
               -outfmt "6 delim=; qseqid sseqid pident"\
               -max_target_seqs 1\
               -qcov_hsp_perc 80
done

# train the amino acid model
TRAIN_FILE=$DATA/train_aa_sintax.fasta
MODEL_FILE=train_aa_$SLURM_ARRAY_TASK_ID

$TIME makeblastdb -in $TRAIN_FILE\
                  -dbtype prot\
                  -out $MODEL_FILE

# classify the amino acid test sequences
for TEST in test testshort
do
TEST_FILE=$DATA/${TEST}_aa.fasta
RAW_FILE=${MODEL}_${TEST}_aa_$SLURM_ARRAY_TASK_ID.raw

$TIME blastp -query $TEST_FILE\
             -db $MODEL_FILE\
             -out $RAW_FILE\
             -num_threads $SLURM_CPUS_PER_TASK\
             -outfmt "6 delim=; qseqid sseqid"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
done

# reformat the output - top hit
for f in *$SLURM_ARRAY_TASK_ID.raw
do
  RESULT_FILE=${f%.raw}.tsv
  RESULT_FILE=$RESULTS/blast_top_hit${f#blast}

  echo "ID	class	Prob_class	order	Prob_order	"\
       "family	Prob_family	subfamily	Prob_subfamily	"\
       "tribe	Prob_tribe	genus	Prob_genus	species	Prob_species"\
       >$RESULT_FILE

  sed -r 's/([^;]+);[^;]+;tax=/\1/; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g'\
      $f\
      >>$RESULT_FILE
done

# reformat the output - thresholds
for f in *_nt_$SLURM_ARRAY_TASK_ID.raw
do
  RESULT_FILE=${f%.raw}.tsv
  RESULT_FILE=$RESULTS/blast_thresh${f#blast}

  echo "ID	class	Prob_class	order	Prob_order	"\
       "family	Prob_family	subfamily	Prob_subfamily	"\
       "tribe	Prob_tribe	genus	Prob_genus	species	Prob_species"\
       >$RESULT_FILE
  awk -F"[,;]"    'BEGIN {
       n=split("90,93,95,96,96.5,97,98", thresh)
     }
     NR>1 {printf "%s", $1
       sub(/tax=/, "")
       for (i = 1; i <= n; i++) {
         if ($10 >= thresh[i]) {
           sub(/[kpcofgs]:/, "", $(i + 2))
           printf "\t%s\t%s", $(i + 2), "1.0"
         } else {
           printf "\tNA\tNA"
         }}
       print ""
     }'\
     $f\
     >>$RESULT_FILE
done

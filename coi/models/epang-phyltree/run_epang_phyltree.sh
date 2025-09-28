#!/usr/bin/env bash
#SBATCH --job-name=iqtree-phyltree
#SBATCH --account=project_2005718
#SBATCH --partition=longrun
#SBATCH --time=14-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --array=1
#SBATCH --output=epang-phyltree_%a.out
#SBATCH --error=epang-phyltree_%a.out
#SBATCH --mail-type=ALL

set -e -x

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
MODEL=epang-phyltree
DATA=../../data
RESULTS=../../results/$MODEL
TRAIN_TAXONOMY_FILE=$DATA/train_gappa.tax

# find IDs of outgroup taxa (arachnids)
awk -F'\t' '/Arachnida/{print $1}' $TRAIN_TAXONOMY_FILE >outgroup.txt

# loop over nt and aa
for ALPHABET in aa nt;
do
  INSTANCE=train_${ALPHABET}
  MODEL_DIR=${INSTANCE}_$SLURM_ARRAY_TASK_ID
  TRAIN_FILE=$DATA/${INSTANCE}_aln.fasta
  CONSTRAINT_FILE=$DATA/phylogenetic_constraints.tre

  rm -fr $MODEL_DIR
  mkdir -p $MODEL_DIR

  if [ "$ALPHABET" = "nt" ]; then
    SUBSTITUTION_MODEL="GTR+G4"
  else
    SUBSTITUTION_MODEL="mtART+G4"
  fi

  # generate reference tree constrained by taxonomic constraint tree
  $TIME iqtree2\
    -s $TRAIN_FILE\
    -T $SLURM_CPUS_PER_TASK\
    -fast\
    -keep-ident\
    --redo\
    -m $SUBSTITUTION_MODEL\
    -g $CONSTRAINT_FILE\
    -pre $MODEL_DIR/$INSTANCE

  # extract substitution model parameters
  if [ "$ALPHABET" = "nt" ]; then
    rates=$(sed -nr '/FINALIZING TREE SEARCH/,$ s|Rate parameters:  A-C: ([0-9.]+) +A-G: ([0-9.]+)  A-T: ([0-9.]+)  C-G: ([0-9.]+)  C-T: ([0-9.]+)  G-T: ([0-9.]+)|\1/\2/\3/\4/\5/\6|p' $MODEL_DIR/$INSTANCE.log)
    echo "Detected rates: $rates"
    frequencies=$(sed -nr '/FINALIZING TREE SEARCH/,$ s|Base frequencies:  A: ([0-9.]+) +C: ([0-9.]+) +G: ([0-9.]+) +T: ([0-9.]+)|\1/\2/\3/\4|p' $MODEL_DIR/$INSTANCE.log)
    echo "Detected base frequencies: $frequencies"
    alpha=$(sed -nr '/FINALIZING TREE SEARCH/,$ s|Gamma shape alpha: +([0-9.]+)|\1|p' $MODEL_DIR/$INSTANCE.log)
    echo "Detected alpha: $alpha"
    SUBSTITUTION_PARAMS="GTR{$rates}+FU{$frequencies}+G4{$alpha}"
  else
    alpha=$(sed -nr '/FINALIZING TREE SEARCH/,$ s|Gamma shape alpha: +([0-9.]+)|\1|p' $MODEL_DIR/$INSTANCE.log)
    echo "Detected alpha: $alpha"
    SUBSTITUTION_PARAMS="mtART+G4{$alpha}"
  fi

  for TEST in test testshort
  do
    TEST_FILE=$DATA/${TEST}_${ALPHABET}_aln.fasta
    RAW_DIR=$MODEL_DIR/$TEST
    mkdir -p $RAW_DIR

    # place test sequences in reference tree
    $TIME epa-ng\
      --tree $MODEL_DIR/$INSTANCE.treefile\
      --ref-msa $TRAIN_FILE\
      --query $TEST_FILE\
      --out-dir $RAW_DIR\
      --model "$SUBSTITUTION_PARAMS"\
      --threads $SLURM_CPUS_PER_TASK\
      --redo

    # assign taxonomy based on phylogenetic placements
    $TIME gappa examine assign\
      --jplace-path $RAW_DIR/epa_result.jplace\
      --taxon-file $TRAIN_TAXONOMY_FILE\
      --root-outgroup outgroup.txt\
      --ranks-string "class|order|family|subfamily|tribe|genus|species"\
      --out-dir $RAW_DIR\
      --file-suffix "_${TEST}_${ALPHABET}_${SLURM_ARRAY_TASK_ID}"\
      --per-query-results\
      --threads $SLURM_CPUS_PER_TASK\
      --allow-file-overwriting
  done
done


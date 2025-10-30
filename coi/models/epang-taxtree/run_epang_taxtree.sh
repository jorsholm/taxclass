#!/usr/bin/env bash
#SBATCH --job-name=epang-taxtree
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
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
MODEL=epang-taxtree
DATA=../../data
RESULTS=../../results/$MODEL
mkdir -p $RESULTS
TRAIN_TAXONOMY_FILE=$DATA/train_gappa.tax

# generate taxonomic constraint tree
$TIME gappa prepare taxonomy-tree\
  --taxon-list-file $TRAIN_TAXONOMY_FILE\
  --out-dir $PWD\
  --threads $SLURM_CPUS_PER_TASK\
  --allow-file-overwriting

# find IDs of outgroup taxa (arachnids)
awk -F'\t' '/Arachnida/{print $1}' $TRAIN_TAXONOMY_FILE >outgroup.txt

# loop over nt and aa
for ALPHABET in aa nt;
do
  INSTANCE=train_${ALPHABET}
  MODEL_DIR=${INSTANCE}_$SLURM_CPUS_PER_TASK
  TRAIN_FILE=$DATA/${INSTANCE}_aln.fasta

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
    -g taxonomy_tree.newick\
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

  for TESTSET in test testshort
  do
    TEST_FILE=$DATA/${TESTSET}_${ALPHABET}_aln.fasta
    RAW_DIR=$MODEL_DIR/$TESTSET
    SUFFIX=${TESTSET}_${ALPHABET}_${SLURM_CPUS_PER_TASK}
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
      --file-suffix "_$SUFFIX"\
      --per-query-results\
      --threads $SLURM_CPUS_PER_TASK\
      --allow-file-overwriting

    RESULT_FILE=$RESULTS/${MODEL}_${SUFFIX}.tsv

    echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	tribe	Prob_tribe	genus	Prob_genus	species	Prob_species" >$RESULT_FILE

    awk -F' ' '
      function write_record(best_prob, best_taxon, new_prob, is_new, new_rank, rank) {
        printf "%s", prev_id
          is_new=0
          new_rank=8
          for (rank = 1; rank <= 7; rank++) {
            if (is_new) {
              printf "\tunk\t%.6f", new_prob[new_rank]
            } else if (new_prob[rank] > best_prob[rank]) {
              is_new = 1
              new_rank = rank
              printf "\tunk\t%.6f", new_prob[rank]
            } else {
              printf "\t%s\t%.6f", best_taxon[rank], best_prob[rank]
            }
          }
          printf "\n"
      }
      BEGIN { prev_id="NA" }
      $1 == "name" { next }
      $1 != prev_id {
        if (prev_id != "NA") {
          write_record(best_prob, best_taxon, new_prob)
        }
        for (rank = 1; rank <= 7; rank++) {
          best_prob[rank] = 0
          new_prob[rank] = 1
          best_taxon[rank] = "NA"
        }
        prev_id = $1
      }
      {
        rank = split($6, a, ";")
        if (rank > 1 && a[rank-1] !~ best_taxon[rank - 1]) next
        new_prob[rank] -= $5
        if ($5 > best_prob[rank]) {
          best_prob[rank] = $5
          best_taxon[rank] = a[rank]
          if (rank < 7) new_prob[rank + 1] = $5
          for (rank2 = rank + 1; rank2 <= 7; rank2++) {
            best_prob[rank2] = 0
            best_taxon[rank2] = "NA"
            new_prob[rank2] = $5
          }
        }
      }
      END{
        write_record(best_prob, best_taxon, new_prob)
      }' "${RAW_DIR}/per_query_${SUFFIX}.tsv" \
      >>$RESULT_FILE

  done
done

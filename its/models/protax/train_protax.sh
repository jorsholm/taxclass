#!/usr/bin/env bash

# Based on readme.txt in protaxA

# 1) set variable NUM_TAXLEVELS based on your taxonomy:
export NUM_TAXLEVELS=7

# 2) add priors for unknown taxa in taxonomy, one value for each level
#    this relates to how much you think there are taxa not included in your taxonomy,
#    larger values add more uncertainty to all predictions
#    unk prior is level-specific where units correspond to leaf nodes of the taxonomy (value * prior(known_species))
#    NOTE: the number of priors in the ,,, list must equal to $NUM_TAXLEVELS

perl $PROTAX/taxonomy_priors.pl "$(cat unk_priors)" taxonomy >taxonomy.priors

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 perl $PROTAX/thintaxonomy.pl $LEVEL taxonomy.priors > tax$LEVEL
done

# 3) generate training data for each level, here 10000 training samples per level

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 echo "LEVEL $LEVEL"
 perl $PROTAX/seqid2taxlevel.pl $LEVEL seqid2tax > ref.tax$LEVEL
 perl $PROTAX/get_all_reference_sequences.pl $LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL
 perl $PROTAX/taxrseq2numeric.pl $LEVEL tax$LEVEL refs.aln > rseqs${LEVEL}.numeric
 perl $PROTAX/generate_training_data.pl $LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL 10000 1 no train$LEVEL 
 perl $PROTAX/traindat2numeric.pl refs.aln train$LEVEL > train${LEVEL}.numeric
done

# to check what kind of training data there is for each level:

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 echo "LEVEL $LEVEL:"
 cut -f4 -d" " train$LEVEL | cut -f1 -d"," | sort | uniq -c
done

# 4) calculate xdat file (sequence similarity predictors), scale the values and save the scaling parameters for later use
#    ...this can take a while if large training data...

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 echo "LEVEL $LEVEL"
 $PROTAX/create_xdata_best2 tax$LEVEL refs.aln rseqs${LEVEL}.numeric train${LEVEL}.numeric > train${LEVEL}.xdat
done

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 perl $PROTAX/scale_xdat.pl sc$LEVEL train${LEVEL}.xdat > train${LEVEL}.scxdat
done

# 5) parameter estimation in R
#    MCMC for parameters separately in each taxonomy level
#    you need to check the convergence and continue iterations or re-initialize adaptive proposal if needed

Rscript ../train_protax.R

#!/usr/bin/env bash

# Based on step1.txt and step2.txt in protax (as well as equivalent code in protaxA)

#############################################
# PART I: Reference sequence database
#############################################

# create USEARCH index:
usearch -makeudb_usearch refs.fa -output refs.udb

#############################################
# PART II: Taxonomy
#############################################


#######################################
# A) Special treatment for Fungi taxonomy: skip nonFungi node in level 1
#######################################

perl $PROTAX/cptaxonomy.pl taxonomy > taxonomy.fungi
perl $PROTAX/thintaxonomy.pl 1 taxonomy.fungi > tax1
# remove nonFungi from taxonomy after kingdom level (1) since it is not divided
# further (this is quick and dirty hack so that no changes need to be done in
# script  generate_training_data.pl)

# note : "dirty hack" and above comment are verbatim from published PROTAX -- BF

#set variable NUM_TAXLEVELS based on your taxonomy:
export NUM_TAXLEVELS=7

for ((LEVEL=2; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 perl $PROTAX/thintaxonomy.pl $LEVEL taxonomy.fungi | grep -v -w nonFungi > tax$LEVEL
done

#############################################
# PART III: Generating training data lists
#############################################

# here 5000 training samples out of which 5% represent unknown taxa
# NOTE: include also LEVEL 1 if normal taxonomy (this is for Fungi example)

for ((LEVEL=2; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 perl $PROTAX/seqid2taxlevel.pl $LEVEL seqid2tax > ref.tax$LEVEL
 perl $PROTAX/get_all_reference_sequences.pl $LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL
 perl $PROTAX/generate_training_data.pl tax$LEVEL ref.tax$LEVEL rseqs$LEVEL 4750 2 no train.level$LEVEL
 perl $PROTAX/generate_unk_training_data.pl $LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL 250 2 no train.unk$LEVEL
 cat train.level$LEVEL train.unk$LEVEL > train$LEVEL
 cut -f6 -d" " train$LEVEL | sort | uniq > train${LEVEL}.id
done

cat train[2,3,4,5,6,7].id | sort | uniq > train.ids

# map training sequences against reference sequences

perl $PROTAX/fastagrep.pl train.ids refs.fa > train.fa
usearch -usearch_global train.fa -db refs.udb -id 0.75 -maxaccepts 1000 -strand both -userfields query+target+id -userout train.m8

for ((LEVEL=2; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 echo $LEVEL
 perl $PROTAX/create_xdata4Q.pl 0.05 train$LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL train.m8 train${LEVEL}.xdat 0
done

# 5) parameter estimation in R
#    MCMC for parameters separately in each taxonomy level
#    you need to check the convergence and continue iterations or re-initialize adaptive proposal if needed

Rscript ../train_protax.R



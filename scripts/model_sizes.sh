#!/usr/bin/env bash

# find sizes of all models
# this is quite idiosyncratic because the models do not store everything in the same file structure

# the model files are not checked in to github, so this needs to be run on a machine which
# has actually run all of the training steps

OUTFILE="results/model_sizes.tsv"
rm -f $OUTFILE
mkdir -p $(dirname $OUTFILE)
touch $OUTFILE

# use extended wildcards
shopt -s extglob

# BayesANT: model is a single rds file
ls -s1 @(coi|its)/models/bayesant/*.rds | sed 's/^ *//' >>$OUTFILE

# BLAST: index is in multiple files, not divided into directories
ls -s1 @(coi|its)/models/blast/*.[np]?? | sed 's/^ *//' >>$OUTFILE

# CREST4: model is in directory
du -s @(coi|its)/models/crest4/train_@(aa|nt)_@(1|4|16|40) | tr "\t" " " >>$OUTFILE

# Dnabarcoder: although various files are produced, only *.best.json is passed for inference
ls -s1 @(coi|its)/models/dnabarcoder/train_@(aa|nt)_@(1|4|16|40)/train_@(aa|nt)_@(1|4|16|40).cutoffs.best.json | sed 's/^ *//' >>$OUTFILE

# IDTAXA: model is a single rds file
ls -s1 @(coi|its)/models/idtaxa/train_@(aa|nt)_@(1|4|16|40).rds | sed 's/^ *//' >>$OUTFILE

# MycoAI: model is a single pt file
ls -s1 @(coi|its)/models/mycoai_@(bert|cnn)/train_nt_@(bert|cnn)_@(1|4|16|40|gpu).pt | sed 's/^ *//' >>$OUTFILE

# Protax-A: model is multiple files in a directory
du -s coi/models/protax-a/train_nt_@(full|train)_tax_1 | tr "\t" " " >>$OUTFILE

# RAXTAX: model is a single bin file
ls -s1 @(coi|its)/models/raxtax/raxtax_@(1|4|16|40)/train_nt_raxtax.bin | sed 's/^ *//' >>$OUTFILE

# RDP-NBC: model is multiple files in a directory
du -s @(coi|its)/models/rdp_nbc/train_nt_1 | tr "\t" " " >>$OUTFILE

# SINTAX: index is a single udb file
ls -s1 @(coi|its)/models/sintax/train_nt_@(1|4|16|40).udb | sed 's/^ *//' >>$OUTFILE

# Protax: model is multiple files in a directory
du -s its/models/protax/train_nt_@(full|train)_tax_@(1|4|16|40) | tr "\t" " " >>$OUTFILE

# EPA-ng: model is two files: the reference tree itself and the IQtree log file
ls -s1 coi/models/epang-@(tax|phyl)tree/train_@(aa|nt)_@(1|4|16|40)/train_@(aa|nt).@(iqtree|treefile) | sed 's/^ *//' >>$OUTFILE



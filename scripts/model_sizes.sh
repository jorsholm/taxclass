#!/usr/bin/env bash

# find sizes of all models
# this is quite idiosyncratic because the models do not store everything in the same file structure

# the model files are not checked in to github, so this needs to be run on a machine which
# has actually run all of the training steps

OUTFILE="results/model_sizes.tsv"
rm -f OUTFILE
mkdir -p $(dirname $OUTFILE)
touch $OUTFILE

# BayesANT: model is a single rds file
ls -s1 */models/bayesant/*.rds | sed 's/^ *//' >>$OUTFILE

# BLAST: index is in multiple files, not divided into directories
ls -s1 */models/blast/*.[np]?? | sed 's/^ *//' >>$OUTFILE

# CREST4: model is in directory
du -s */models/crest4/train_??_* | tr "\t" " " >>$OUTFILE

# Dnabarcoder: although various files are produced, only *.best.json is passed for inference
ls -s1 */models/dnabarcoder/train_??_*/train_??_*.cutoffs.best.json | sed 's/^ *//' >>$OUTFILE

# IDTAXA: model is a single rds file
ls -s1 */models/idtaxa/train_??_.*.rds | sed 's/^ *//' >>$OUTFILE

# MycoAI: model is a single pt file
ls -s1 */models/mycoai_*/train_??_*.pt | sed 's/^ *//' >>$OUTFILE

# Protax-A: model is multiple files in a directory
du -s */models/protax-a/train_nt_*_tax_1 | tr "\t" " " >>$OUTFILE

# RDP-NBC: model is multiple files in a directory
du -s */models/rdp_nbc/train_nt_1 | tr "\t" " " >>$OUTFILE

# SINTAX: index is a single udb file
ls -s1 */models/sintax/train_??_*.udb | sed 's/^ *//' >>$OUTFILE

# ProtaxFungi: TBD

# EPA-ng: TBD



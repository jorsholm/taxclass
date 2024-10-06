#!/bin/sh

SRC="97_otus"
SEQFILE="${SRC}.fasta"
TREEFILE="${SRC}.tree"
TAXONFILE="${SRC}_taxonomy.txt"

MODEL="GTR"
DBNAME="gg_97_otus"
DB="${DBNAME}_${MODEL}"

echo "Buiding database $DB"

hmmufotu-build $SEQFILE $TREEFILE -a $TAXONFILE -n $DBNAME -s $MODEL -p 4 -v

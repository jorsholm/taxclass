#!/usr/bin/env bash

for f in */models/*/*.out
do
 amplicon=$(dirname $(dirname $(dirname $f)))
 model=$(basename $(dirname $f))
 awk -f scripts/extract_stats.awk <$f >$amplicon/results/$model/$(basename $f .out)_time.tsv
done

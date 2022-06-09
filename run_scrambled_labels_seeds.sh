#!/bin/bash

outdir="scrambled_labels_seeds"

mkdir -p data/$outdir

for seed in {-1..9}
do
  echo "Preparing PNET input data ..."
  Rscript prepare_scrambled_labels.R FALSE $seed

  echo "Running PNET ..."
  python main.py -1 $outdir

  mv data/$outdir/234_20080808 data/$outdir/$seed
done

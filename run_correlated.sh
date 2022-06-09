#!/bin/bash

outdir="correlated"

mkdir -p data/$outdir

echo "Preparing PNET input data ..."
Rscript prepare_correlated.R

echo "Running PNET ..."
for seed in {-1..9}
do
  python main.py $seed $outdir
done

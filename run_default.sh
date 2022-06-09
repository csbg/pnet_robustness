#!/bin/bash

outdir="default"

mkdir -p data/$outdir

echo "Preparing PNET input data ..."
Rscript prepare_default.R

echo "Running PNET ..."
for seed in {-1..9}
do
  python main.py $seed $outdir
done

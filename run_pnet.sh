#!/bin/bash

experiment_name=${1:-default}

results_dir="data/$experiment_name"
mkdir -p $results_dir

echo "Running PNET ..."
for seed in {-1..49}
do
  if [ $seed -eq -1 ]
  then
    random_seed="234"
    pipeline_seed="20080808"
  else
    random_seed=$seed
    pipeline_seed=$seed
  fi

  python pnet_prostate_paper/train/run_me.py $random_seed $pipeline_seed
  python pnet_prostate_paper/analysis/prepare_data.py

  target_dir="$results_dir/${random_seed}_${pipeline_seed}"
  mkdir -p $target_dir/onsplit_average_reg_10_tanh_large_testing/fs
  cp pnet_prostate_paper/analysis/extracted/* $target_dir
  cp pnet_prostate_paper/_logs/p1000/pnet/onsplit_average_reg_10_tanh_large_testing/fs/* $target_dir/onsplit_average_reg_10_tanh_large_testing/fs/
  for file in all_scores.csv log.log P-net_ALL_params.yml P-net_ALL_testing.csv P-net_ALL_training.csv
  do
    cp pnet_prostate_paper/_logs/p1000/pnet/onsplit_average_reg_10_tanh_large_testing/$file $target_dir/onsplit_average_reg_10_tanh_large_testing/
  done
done

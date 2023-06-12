#!/bin/bash

# convenience logging function
log() {
  echo "pnet_robustness [$(date -Is)]:" "$@"
}

# setup environment
source /opt/conda/etc/profile.d/conda.sh
conda activate pnet_env
export PYTHONPATH=/app/pnet_prostate_paper:$PYTHONPATH

# create directory for results
pnet_dir=/app/pnet_prostate_paper/
results_dir="data/$1"
mkdir -p $results_dir

# run PNET with several seeds
log "Start PNET analysis"
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

  log "Execute run_me.py"
  python $pnet_dir/train/run_me.py $random_seed $pipeline_seed

  log "Execute prepare_data.py"
  python $pnet_dir/analysis/prepare_data.py

  # copy results
  log "Save results"
  target_dir="$results_dir/${random_seed}_${pipeline_seed}"
  mkdir -p $target_dir

  cp $pnet_dir/analysis/extracted/node_importance_graph_adjusted.csv $target_dir/node_importance.csv
  cp $pnet_dir/_logs/p1000/pnet/onsplit_average_reg_10_tanh_large_testing/P-net_ALL_testing.csv $target_dir/predictions_test.csv
  cp $pnet_dir/_logs/p1000/pnet/onsplit_average_reg_10_tanh_large_testing/P-net_ALL_training.csv $target_dir/predictions_train.csv

  # during the initial run, also save files describing the neural net structure
  if [ $seed -eq -1 ]
  then
    cp $pnet_dir/analysis/extracted/link_weights_?.csv $target_dir/
  fi
done

log "Done"

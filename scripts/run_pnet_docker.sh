#!/bin/bash

docker run \
  --rm -it \
  -v $PWD:/pnet_robustness \
  -v $PWD/pnet_data/mounted/mutations.csv:/app/pnet_prostate_paper/_database/prostate/processed/P1000_final_analysis_set_cross_important_only.csv \
  -v $PWD/pnet_data/mounted/cnvs.csv:/app/pnet_prostate_paper/_database/prostate/processed/P1000_data_CNA_paper.csv \
  -v $PWD/pnet_data/mounted/labels.csv:/app/pnet_prostate_paper/_database/prostate/processed/response_paper.csv \
  -v $PWD/pnet_data/mounted/training_set.csv:/app/pnet_prostate_paper/_database/prostate/splits/training_set_0.csv \
  -v $PWD/pnet_data/mounted/test_set.csv:/app/pnet_prostate_paper/_database/prostate/splits/test_set.csv \
  -v $PWD/pnet_data/mounted/validation_set.csv:/app/pnet_prostate_paper/_database/prostate/splits/validation_set.csv \
  ghcr.io/csbg/pnet-container:1.0.0 $@

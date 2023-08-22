#!/bin/bash

singularity exec \
  --writable-tmpfs \
  --bind $PWD:/pnet_robustness \
  --bind $PWD/pnet_data/mounted/mutations.csv:/app/pnet_prostate_paper/_database/prostate/processed/P1000_final_analysis_set_cross_important_only.csv \
  --bind $PWD/pnet_data/mounted/cnvs.csv:/app/pnet_prostate_paper/_database/prostate/processed/P1000_data_CNA_paper.csv \
  --bind $PWD/pnet_data/mounted/labels.csv:/app/pnet_prostate_paper/_database/prostate/processed/response_paper.csv \
  --bind $PWD/pnet_data/mounted/training_set.csv:/app/pnet_prostate_paper/_database/prostate/splits/training_set_0.csv \
  --bind $PWD/pnet_data/mounted/test_set.csv:/app/pnet_prostate_paper/_database/prostate/splits/test_set.csv \
  --bind $PWD/pnet_data/mounted/validation_set.csv:/app/pnet_prostate_paper/_database/prostate/splits/validation_set.csv \
  pnet-container_1.0.0.sif /app/entrypoint.sh $@

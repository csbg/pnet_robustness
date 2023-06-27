#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate dtox_env
export PYTHONPATH=/app/DTox/code:$PYTHONPATH
python /app/run_dtox.py $@

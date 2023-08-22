#!/bin/bash

singularity exec \
  --writable-tmpfs \
  --bind $PWD:/dtox_robustness \
  dtox-container_1.0.0.sif /app/entrypoint.sh $@ $@

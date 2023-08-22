#!/bin/bash

singularity run \
  --writable-tmpfs \
  --bind $PWD:/pnet_robustness \
  r_pnet_robustness_1.0.0.sif $@

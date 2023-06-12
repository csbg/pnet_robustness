#!/bin/bash

docker run \
  --rm -it \
  -v $PWD:/pnet_robustness \
 ghcr.io/csbg/pnet-container $@

# docker run --rm -it -v $PWD:/pnet_robustness ghcr.io/csbg/pnet-container ./run_pnet.sh original

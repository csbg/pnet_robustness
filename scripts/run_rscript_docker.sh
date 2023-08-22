#!/bin/bash

docker run \
  --rm -it \
  -v $PWD:/pnet_robustness \
  ghcr.io/csbg/r_pnet_robustness:1.0.0 $@

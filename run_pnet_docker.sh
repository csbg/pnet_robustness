#!/bin/bash

docker run \
  --rm -it \
  -v $PWD:/pnet_robustness \
 ghcr.io/csbg/pnet-container $@

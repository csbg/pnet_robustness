#!/bin/bash

docker run \
  --rm -it \
  -v $PWD:/dtox_robustness \
  ghcr.io/csbg/dtox-container:1.0.0 $@

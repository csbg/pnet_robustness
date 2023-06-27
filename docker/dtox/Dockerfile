# syntax=docker/dockerfile:1

FROM continuumio/miniconda3:23.3.1-0

WORKDIR /app
COPY entrypoint.sh environment_dtox.yml patch_seeds.diff run_dtox.py setup.sh ./
RUN bash setup.sh

WORKDIR /dtox_robustness
SHELL ["/bin/bash", "--login", "-c"]
ENTRYPOINT ["/app/entrypoint.sh"]

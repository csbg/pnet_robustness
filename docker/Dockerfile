# syntax=docker/dockerfile:1

FROM continuumio/miniconda3:23.3.1-0

WORKDIR /app
COPY entrypoint.sh environment_pnet.yml patch_seeds.diff setup.sh ./
RUN bash setup.sh

WORKDIR /pnet_robustness
ENV MPLBACKEND=Agg
SHELL ["/bin/bash", "--login", "-c"]
ENTRYPOINT ["/app/entrypoint.sh"]

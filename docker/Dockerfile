# docker build --tag ghcr.io/csbg/pnet-container .

FROM continuumio/miniconda3:23.3.1-0

WORKDIR /app
COPY environment_pnet.yml patch_seeds.diff run_pnet.sh setup.sh ./
RUN bash setup.sh

WORKDIR /pnet_robustness
ENV MPLBACKEND=Agg
SHELL ["/bin/bash", "--login", "-c"]
ENTRYPOINT ["/app/run_pnet.sh"]

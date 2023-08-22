# syntax=docker/dockerfile:1

FROM rocker/r-ver:4.3.1
WORKDIR /app
COPY renv.lock ./
RUN apt-get update
RUN apt-get install -y --no-install-recommends libglpk-dev
RUN R -e 'install.packages("renv"); \
          renv::consent(provided = TRUE); \
          renv::settings$use.cache(FALSE); \
          renv::restore(lockfile = "/app/renv.lock", \
                        repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))'
WORKDIR /pnet_robustness
ENTRYPOINT ["Rscript", "--vanilla"]

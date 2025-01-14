# Reliable interpretability of biology-inspired deep neural networks


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8386694.svg)](https://doi.org/10.5281/zenodo.8386694)

This code supplements the [publication](https://doi.org/10.1038/s41540-023-00310-8) by Esser-Skala and Fortelny (2023).


## Folders

(Not all of these folders are included in the git repository.)

- `data`: output files generated by P-NET and DTox (available from https://doi.org/10.5281/zenodo.7760561)
- `doc`: project documentation
- `docker`: files for creating [Docker containers](#appendix-how-to-build-docker-images) with P-NET or DTox installed
- `literature`: relevant publications
- `plots`: generated plots
- `pnet_data`: P-NET data files
- `renv`: R environment data
- `scripts`: bash and R scripts



## Preparation

### P-NET and DTox

Download the provided Docker containers from the GitHub Container registry:

```bash
docker pull ghcr.io/csbg/pnet-container:1.0.0
docker pull ghcr.io/csbg/dtox-container:1.0.0
```

Alternatively, pull these containers with Apptainer/Singularity:

```bash
singularity pull docker://ghcr.io/csbg/pnet-container:1.0.0
singularity pull docker://ghcr.io/csbg/dtox-container:1.0.0
```

When using the latter container format, replace all calls to `run_[pnet/dtox]_docker.sh` with `run_[pnet/dtox]_singularity.sh`.


### R

In order to run the R scripts, you will need one of the following:

- an installation of R 4.3.1; restore required packages from `renv.lock` via
  ```bash
  Rscript -e "renv::restore()"`
  ```

- the Docker container available from the GitHub Container registry:
  ```bash
  docker pull ghcr.io/csbg/r_pnet_robustness:1.0.0
  ```
  Replace calls to `Rscript` below by `scripts/run_rscript_docker.sh`.

- the Apptainer/Singularity container:
  ```bash
  singularity pull docker://ghcr.io/csbg/r_pnet_robustness:1.0.0
  ```
  Replace calls to `Rscript` below by `scripts/run_rscript_singularity.sh`.



### Datasets

Download the [MSK-IMPACT 2017](https://www.nature.com/articles/nm.4333) dataset:

```bash
wget https://cbioportal-datahub.s3.amazonaws.com/msk_impact_2017.tar.gz
tar xzf msk_impact_2017.tar.gz -C pnet_data
```



## Run P-NET experiments

Generally, each experiment comprises the following steps:

1. Load P-NET input data via `load_data_[dataset].R`.
2. Optionally, modify input data via `modify_data_[technique].R`.
3. Run P-NET via Docker using the provided bash script `run_pnet_docker.sh`. This script has three arguments:
  - `-e experiment`: experiment name, required
  - `-l [n]`: lower seed, optional (default: -1, which uses the original seeds)
  - `-u [n]`: upper seed, optional (default: 49)

Within each experiment, results from each run are saved in a subfolder indicating the two random seeds used (e.g., `data/pnet_original/0_0`).

`utils.R` is required by all data preparation scripts.


### Original setup

Run P-NET with the original setup as described in the publication.

```bash
Rscript scripts/load_data_original.R
scripts/run_pnet_docker.sh -e pnet_original
```


### Deterministic inputs

Input data is modified so that presence of mutation and copy number amplification is perfectly correlated with class label 1 (copy number deletion is always 0).

```bash
Rscript scripts/load_data_original.R
Rscript scripts/modify_data_deterministic.R
scripts/run_pnet_docker.sh -e pnet_deterministic
```


### Shuffled labels

Shuffle training/test labels before each run using uniform class frequencies.

```bash
for seed in {-1..49}; do
  Rscript scripts/load_data_original.R
  Rscript scripts/modify_data_shuffled.R FALSE $seed
  scripts/run_pnet_docker.sh -e pnet_shuffled_each -l $seed -u $seed
done
```


### MSK-IMPACT 2017 dataset

```bash
Rscript scripts/load_data_mskimpact.R "Non-Small Cell Lung Cancer"
scripts/run_pnet_docker.sh -e mskimpact_nsclc_original

for seed in {-1..49}; do
  Rscript scripts/load_data_mskimpact.R "Non-Small Cell Lung Cancer"
  Rscript scripts/modify_data_shuffled.R FALSE $seed
  scripts/run_pnet_docker.sh -e mskimpact_nsclc_shuffled -l $seed -u $seed
done


Rscript scripts/load_data_mskimpact.R "Breast Cancer"
scripts/run_pnet_docker.sh -e mskimpact_bc_original

for seed in {-1..49}; do
  Rscript scripts/load_data_mskimpact.R "Breast Cancer"
  Rscript scripts/modify_data_shuffled.R FALSE $seed
  scripts/run_pnet_docker.sh -e mskimpact_bc_shuffled -l $seed -u $seed
done


Rscript scripts/load_data_mskimpact.R "Colorectal Cancer"
scripts/run_pnet_docker.sh -e mskimpact_cc_original

for seed in {-1..49}; do
  Rscript scripts/load_data_mskimpact.R "Colorectal Cancer"
  Rscript scripts/modify_data_shuffled.R FALSE $seed
  scripts/run_pnet_docker.sh -e mskimpact_cc_shuffled -l $seed -u $seed
done


Rscript scripts/load_data_mskimpact.R "Prostate Cancer"
scripts/run_pnet_docker.sh -e mskimpact_pc_original

for seed in {-1..49}; do
  Rscript scripts/load_data_mskimpact.R "Prostate Cancer"
  Rscript scripts/modify_data_shuffled.R FALSE $seed
  scripts/run_pnet_docker.sh -e mskimpact_pc_shuffled -l $seed -u $seed
done
```



## Run DTox experiments

Run DTox with seeds ranging from 0 (i.e., the original seed) to 50.

```bash
scripts/run_dtox_docker.sh
```

Results from each run are saved in a subfolder indicating the random seed used (e.g., `data/dtox/0`).



## Analyze results

`plot_figures.R` generates all figures shown in the publication, using files in `data` (described below):

```bash
Rscript scripts/plot_figures.R
```

`styling.R` is required by this script.


### P-NET

After each run, the following files are copied from the P-NET output folders:

- `analysis/extracted/node_importance_graph_adjusted.csv` (renamed to `node_importance.csv`): contains node importance scores, with the following columns:
  - (first, unnamed): node name
  - coef: original node importance scores
  - coef_graph: indegree plus outdegree of node
  - coef_combined: adjusted node importance score (= coef / coef_graph if coef_graph > mean(coef_graph) + 5 sd(coef_graph) in the respective layer)
  - coef_combined_zscore: scaled coef_combined
  - coef_combined2: z(z(coef_graph) - z(coef))
  - layer: layer of the node
- `_logs/p1000/pnet/onsplit_average_reg_10_tanh_large_testing/P-net_ALL_testing.csv` (renamed to `predictions_test.csv`): predictions for the test set, with the following columns:
  - (first, unnamed): sample name
  - pred: predicted class (unfortunately, encoded by a double 1.0 or 0.0)
  - pred_scores: probability of the predicted class
  - y: true class (encoded as integer 1 or 0)
- `_logs/p1000/pnet/onsplit_average_reg_10_tanh_large_testing/P-net_ALL_training.csv` (renamed to `predictions_train.csv`): predictions for the training set (same columns as above)


### DTox

The following files generated by DTox are required for subsequent analyses:

- `module_relevance.tsv`: contains node importance scores, with the following columns:
  - (first, unnamed): compound identifier
  - remaining columns: node identifiers (UniProt and Reactome IDs)
- `test_labels.csv`: predictions for the test set, with two columns:
  - truth: true label (0 or 1)
  - predicted: predicted label (decimal number between 0 and 1)



## Appendix: How to build Docker images

### P-NET

The folder `docker/pnet` contains everything needed for building a Docker image with P-NET installed:

- `Dockerfile`: instructions for assembling the image
- `entrypoint.sh`: script for running P-NET; used as entrypoint in the container
- `environment_pnet.yml`: conda environment specification
- `patch_seeds.diff`: patch that allows to change the random seed for P-NET
- `setup.sh`: executed during image assembly; installs P-NET (with input data) and conda


Build and deploy this image via

```bash
docker build --tag ghcr.io/csbg/pnet-container:1.0.0 .
docker push ghcr.io/csbg/pnet-container:1.0.0
```


### DTox

The folder `docker/dtox` contains everything needed for building a Docker image with DTox installed:

- `Dockerfile`: instructions for assembling the image
- `entrypoint.sh`: entrypoint in the container; activates conda environment
- `environment_dtox.yml`: conda environment specification
- `patch_seeds.diff`: patch that allows to change the random seed for DTox and saves predicted labels for the test set
- `run_dtox.py`: executes the DTox workflow as described in the [tutorial](https://github.com/EpistasisLab/DTox/blob/master/DTox-implementation.ipynb) available in the DTox GitHub repository
- `setup.sh`: executed during image assembly; installs DTox and conda


Build and deploy this image via

```bash
docker build --tag ghcr.io/csbg/dtox-container:1.0.0 .
docker push ghcr.io/csbg/dtox-container:1.0.0
```


### R

The folder `docker/r` contains the `Dockerfile` needed for building a Docker image with R and required packages installed.

Build and deploy this image via

```bash
cp ../../renv.lock .
docker build --tag ghcr.io/csbg/r_pnet_robustness:1.0.0 .
docker push ghcr.io/csbg/r_pnet_robustness:1.0.0
```



## Appendix: Description of files required by P-NET

- `genes/`: only the genes present in both of the following two files will be analyzed:
            (a) `tcga_prostate_expressed_genes_and_cancer_genes.csv`
            (b) `HUGO_genes/protein-coding_gene_with_coordinate_minimal.txt`
                 (TSV, no column names; meaning of columns: chromosome, start, end, gene name)
- `pathways/`:
  - `pathways_short_names.xlsx`: short pathway names for figure labels
  - `Reactome/ReactomePathways.gmt`:
    genes associated with Reactome pathways;
    TSV, no column names, variable number of columns:
    (1) pathway name
    (2) reactome id
    (3) type (unused)
    (4ff) associated genes
  - `Reactome/ReactomePathways.txt`:
    TSV mapping Reactome ids to names;
    loaded by P-NET but apparently not used (?)
  - `Reactome/ReactomePathwaysRelation.txt`:
    TSV specifying the Reactome pathway hierarchy as edge list;
    columns indicate parent and child;
    only human pathways are used (i.e., the child id has to start with "HSA")
- `prostate/`
  - `processed/`
    - `P1000_final_analysis_set_cross_important_only.csv`: mutation data;
      first column contains sample name,
      remaining columns represent genes,
      cells indicate number of mutations;
      data is preprocessed to a binary matrix, indicating presence/absence
      of at least one mutation (i.e., 1 if original >= 1)
    - `P1000_data_CNA_paper.csv`: CNV data;
      first column (unnamed) contains sample name,
      remaining columns represent genes,
      cells indicate copy number status;
      data is preprocessed to two binary matrices:
      one indicates presence of copy number amplification (1 if original > 1.5),
      the other indicates presence of CN deletion (1 if original < -1.5)
    - `response_paper.csv`: input labels, two columns:
                            (1) id – sample name
                            (2) response – sample label (1 = metastatic tumor)
  - `splits/`: splits of input data; all files have three columns:
               (1) [unnamed] – running number starting at zero
               (2) id – sample name
               (3) response – sample label (column is NOT used by P-NET!)
    - `test_set.csv`: samples in the test set
    - `training_set_0.csv`: samples in the training set
    - `validation_set.csv`: samples in the validation set

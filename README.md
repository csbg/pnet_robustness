# PNET robustness

## Prepare PNET

Clone the git repository and change `pipeline/one_split.py` and `train/run_me.py` to allow selection of random seeds.

```bash
gh repo clone marakeby/pnet_prostate_paper
cd pnet_prostate_paper
git checkout -b pnet-robustness 2b16264
git apply ../pnet_data/patch_seeds.diff
git commit -a -m "changes code for PNET robustness tests"
conda env create --name pnet_env --file=environment.yml
```

Download [data files](https://drive.google.com/uc?id=17nssbdUylkyQY1ebtxsIw5UzTAd0zxWb&export=download) and move the decompressed folder `_database` into `pnet_prostate_paper`. PNET requires the following files:

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
    loaded by PNET but apparently not used (?)
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
               (3) response – sample label (column is NOT used by PNET!)
    - `test_set.csv`: samples in the test set
    - `training_set_0.csv`: samples in the training set
    - `validation_set.csv`: samples in the validation set



## Run experiments

Set up the environment:

``` bash
# first line only required if conda is not activated in .bashrc
source /home/wolfgang/Programs/miniconda3/etc/profile.d/conda.sh
conda activate pnet_env
export PYTHONPATH=$PWD/pnet_prostate_paper:$PYTHONPATH
```

Generally, each experiment comprises two steps:
1. Prepare PNET input data via `prepare_*.R`.
2. Run PNET via `run_pnet.sh [experiment name]` (the argument may be empty, in which case it has the value `default`).

Within each experiment, results from each run are saved in a subfolder indicating the two random seeds used (e.g., `0_0`).


### Default settings

Run PNET with default settings.

```bash
Rscript prepare_default.R
./run_pnet.sh default
```


### Correlated predictors and classes

Input data is modified so that presence of mutation and copy number amplification is perfectly correlated with class label 1 (copy number deletion is always 0).

```bash
Rscript prepare_correlated.R
./run_pnet.sh correlated
```


### Dropout

(1) Set dropout to 0.95 by modifying lines 21 and 37 in `train/params/P1000/pnet/onsplit_average_reg_10_tanh_large_testing.py` (revert changes afterwards).

```bash
git apply --directory=pnet_prostate_paper pnet_data/patch_high_dropout.diff
./run_pnet.sh dropout_high
git apply -R --directory=pnet_prostate_paper pnet_data/patch_high_dropout.diff
```

(2) Eliminate dropout layers altogether by removing lines 167, 169, 236, and 237 in `model/builders/builders_utils.py` (revert changes afterwards).

```bash
git apply --directory=pnet_prostate_paper pnet_data/patch_no_dropout.diff
./run_pnet.sh dropout_none
git apply -R --directory=pnet_prostate_paper pnet_data/patch_no_dropout.diff
```


### Scrambled labels

(1) Scramble training/test labels while retaining class frequency.

```bash
Rscript prepare_scrambled_labels.R TRUE 0
./run_pnet.sh scrambled_labels_original_class_frequency
```

(2) Use uniform class frequencies.

```bash
Rscript prepare_scrambled_labels.R FALSE 0
./run_pnet.sh scrambled
```


### Scrambled features

Scramble input features. Values in the mutation (CNA) matrix are replaced by 0 and 1 (2); the nonzero element is chosen with a given probability (50%, 5%, or 0.1%). A random seed of 0 or 1 is used for scrambling.

```bash
Rscript prepare_scrambled_features.R 0.5 0
./run_pnet.sh scrambled_features_0.5_seed_0

Rscript prepare_scrambled_features.R 0.05 0
./run_pnet.sh scrambled_features_0.05_seed_0

Rscript prepare_scrambled_features.R 0.001 0
./run_pnet.sh scrambled_features_0.001_seed_0

Rscript prepare_scrambled_features.R 0.5 1
./run_pnet.sh scrambled_features_0.5_seed_1

Rscript prepare_scrambled_features.R 0.05 1
./run_pnet.sh scrambled_features_0.05_seed_1
```



## PNET output data

These files are directly copied from the PNET output folders (`analysis/extracted` and `_logs/p1000/pnet`):

- `node_importance_graph_adjusted.csv` contains node importance scores; important columns:
  - (first, unnamed): node name
  - coef: original node importance scores
  - coef_graph: indegree plus outdegree of node
  - coef_combined: adjusted node importance score (= coef / coef_graph if coef_graph > mean(coef_graph) + 5 sd(coef_graph) in the respective layer)
  - coef_combined_zscore: scaled coef_combined
  - coef_combined2: z(z(coef_graph) - z(coef))
  - layer: layer of the node
- `link_weights_X.csv`: link weights from layer X (nodes in row) to the next layer (nodes as columns)
- `onsplit_average_reg_10_tanh_large_testing/P-net_ALL_testing.csv`: predictions for the test set, with the following columns:
  - (first, unnamed): sample name
  - pred: predicted class (unfortunately, encoded a double 1.0 or 0.0)
  - pred_scores: probability of the predicted class
  - y: true class (encoded as integer 1 or 0)



## Sync with GFS

```bash
rsync -azvhP --delete /home/wolfgang/Plus/Projects/pnet_robustness/data/ /mnt/agfortelny/people/wskala/pnet_robustness/data
```

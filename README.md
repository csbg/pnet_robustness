# PNET robustness

## Prepare PNET

Clone the git repository and change `pipeline/one_split.py` and `train/run_me.py` to allow selection of random seeds (use patch in `patch/`).

```bash
gh repo clone marakeby/pnet_prostate_paper
git checkout 2b16264
git apply pnet_patch.diff
```

Download data files from [https://drive.google.com/uc?id=17nssbdUylkyQY1ebtxsIw5UzTAd0zxWb&export=download] and unzip the folder `_database`. (Note: Only files described [below](#pnet-input-data) are required.)



## Run an experiment

Execute all of the following commands in the root folder (i.e., `pnet_prostate_cancer`).

Set up the environment via

``` bash
eval /home/wolfgang/Programs/miniconda3/bin/conda "shell.fish" "hook" $argv | source
conda activate pnet_env
set -x PYTHONPATH $PWD
```

Run the respective experiment via one of the `seed_tests/run_*.sh` scripts, e.g.,

``` bash
seed_tests/run_default.sh
```

or, if only the default seed is required, via

```bash
python seed_tests/main.py
```


## PNET input data

PNET employs the following data files, all of which are stored in `_database/`:

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


The three CSV files in `prostate/processed/` have been renamed to `*.csv.original` with read-only permission.



## Experiments

Prior to each experiment, `run_*.sh` executes the respective `prepare_*.R` script to prepare PNET input data.

Input data for this R script is saved in subfolders of `data/`:

- `default`:
  run PNET with default settings (`run_default.sh`)
- `correlated`:
  input data is modified so that presence of mutation and copy number
  amplification is perfectly correlated with class label 1
  (`run_correlated.sh`); copy number deletion is always 0
- `dropout_high`:
  set dropout to 0.95 by modifying lines 21 and 37 in
  `train/params/P1000/pnet/onsplit_average_reg_10_tanh_large_testing.py`
  (run via `run_default.R`; do not forget to revert these changes!)
- `dropout_none`:
  eliminate dropout layers (i.e., comment out lines 167, 169 and 236,
  and 237 in `model/builders/builders_utils.py`)
  (run via `run_default.R`; do not forget to revert these changes!)
- `scrambled_labels`: scramble training/test labels
  while retaining class frequency (`run_scrambled_labels.R`)
- `scrambled_labels_balanced`:
  as above, but equal class frequencies (`run_scrambled_labels_balanced.R`)
- `scrambled_features_0.5_seed_0`
- `scrambled_features_0.05_seed_0`
- `scrambled_features_0.001_seed_0`
- `scrambled_features_0.5_seed_1`
- `scrambled_features_0.05_seed_1`:
  scramble input features (`run_scrambled_features_*_seed_*.R`);
  values in the mutation (CNA) matrix are replaced by 0 and 1 (2);
  the nonzero element is chosen with the given probability (50%, 5%, or 0.1%);
  a random seed of 0 or 1 is used from scrambling



## PNET output data

These files are directly copied from the PNET output folders (`analysis/extracted` and `_logs/p1000/pnet`).

- Results from each run are saved in a subfolder indicating the two random seeds used (e.g., `0_0`).
- `node_importance_graph_adjusted.csv` contains node importance scores; important columns:
  - (first, unnamed): node name
  - coef: node importance scores
  - layer: layer of the node
- `link_weights_X.csv`: link weights from layer X (nodes in row) to the next layer (nodes as columns)
- `onsplit_average_reg_10_tanh_large_testing/P-net_ALL_testing.csv`: predictions for the test set, with the following columns:
  - (first, unnamed): sample name
  - pred: predicted class (unfortunately, encoded a double 1.0 or 0.0)
  - pred_scores: probability of the predicted class
  - y: true class (encoded as integer 1 or 0)


You may remove unnecessary files via

```bash
# cd data/[experiment]
rm -r */*/*.png */*/Logistic* */*/pnet*history.csv */*/pnet*validation
```


## Sync with GFS

```bash
rsync -azvhP --delete /home/wolfgang/Plus/Experiments/pnet_prostate_paper/seed_tests/data/ /mnt/agfortelny/people/wskala/pnet_prostate_paper/seed_tests/data
```

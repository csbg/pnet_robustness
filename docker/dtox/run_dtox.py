"""Run DTox with different seeds"""

import logging
import os
from sys import argv
import pandas as pd
import dtox
import dtox_interpret
import targettox

logging.basicConfig(format="%(levelname)s [%(asctime)s] %(message)s",
                    level=logging.INFO)


# constants
AO_COL = 'assay_outcome'
START_DIR = os.getcwd()
os.chdir("/app/DTox")

if len(argv) != 3:
    LOWER_SEED = 0
    UPPER_SEED = 50
else:
    LOWER_SEED = int(argv[1])
    UPPER_SEED = int(argv[2])
logging.info("Using seeds %i to %i", LOWER_SEED, UPPER_SEED)

# (1) load raw input data
logging.info("Loading input data")
train_data_raw = pd.read_csv(
    'data/example/mitotox_example_input_train.tsv',
    sep='\t',
    header=0,
    index_col=0
)
test_data_raw = pd.read_csv(
    'data/example/mitotox_example_input_test.tsv',
    sep='\t',
    header=0,
    index_col=0
)

# (2) infer binding profiles (actual training/test data)
logging.info("Infer target binding profile")
train_data = targettox.derive_target_profile(train_data_raw)
train_data[AO_COL] = train_data_raw[AO_COL]
test_data = targettox.derive_target_profile(test_data_raw)
test_data[AO_COL] = test_data_raw[AO_COL]

# (3) train several models and get importance scores
for seed in range(LOWER_SEED, UPPER_SEED + 1):
    result_folder = f"{START_DIR}/data/dtox/{seed}/"
    os.makedirs(result_folder, exist_ok=True)

    logging.info("Train model with seed %s", seed)
    model_info, model, model_loss, model_training_summary = dtox.dtox(
        train_data,
        label_col_name=AO_COL,
        out_folder=result_folder,
        seed=seed
    )

    logging.info("Evaluate model performance")
    train_performance, test_performance = dtox.dtox_eval(
        train_data,
        test_data,
        label_col_name=AO_COL,
        trained_model=model,
        loss=model_loss,
        out_folder=result_folder
    )

    logging.info("Calculate node importances")
    pathway_rel_df, rel_fdr_df = dtox_interpret.dtox_interpret(
        test_data[test_data[AO_COL] == 1],
        dtox_model=model,
        dtox_combine_df=train_data,
        label_col_name=AO_COL,
        out_folder=result_folder,
        N_null_models=1
    )

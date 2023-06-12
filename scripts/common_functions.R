# common functions and variables

MOUNTED_FILES <- list(
  mutations = "pnet_data/mutations.csv",
  cnvs = "pnet_data/cnvs.csv",
  labels = "pnet_data/labels.csv",
  training_set = "pnet_data/training_set.csv",
  test_set = "pnet_data/test_set.csv",
  validation_set = "pnet_data/validation_set.csv"
)

ORIGINAL_FILES <- list(
  mutations = "pnet_data/original/P1000_final_analysis_set_cross_important_only.csv",
  cnvs = "pnet_data/original/P1000_data_CNA_paper.csv",
  labels = "pnet_data/original/response_paper.csv",
  training_set = "pnet_data/original/training_set_0.csv",
  test_set = "pnet_data/original/test_set.csv",
  validation_set = "pnet_data/original/validation_set.csv"
)

restore_input_files <- function() {
  walk2(
    ORIGINAL_FILES,
    MOUNTED_FILES,
    file_copy,
    overwrite = TRUE
  )
}

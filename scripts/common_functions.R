# common functions and variables

MOUNTED_FILES <- list(
  mutations = "pnet_data/mounted/mutations.csv",
  cnvs = "pnet_data/mounted/cnvs.csv",
  labels = "pnet_data/mounted/labels.csv",
  training_set = "pnet_data/mounted/training_set.csv",
  test_set = "pnet_data/mounted/test_set.csv",
  validation_set = "pnet_data/mounted/validation_set.csv"
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
  dir_create("pnet_data/mounted")
  walk2(
    ORIGINAL_FILES,
    MOUNTED_FILES,
    file_copy,
    overwrite = TRUE
  )
}

# common functions and variables

MUT_FILE <- "../_database/prostate/processed/P1000_final_analysis_set_cross_important_only.csv.original"
CNV_FILE <- "../_database/prostate/processed/P1000_data_CNA_paper.csv.original"
LABEL_FILE <- "../_database/prostate/processed/response_paper.csv.original"

restore_input_file <- function(original_file) {
  copied_file <- path_ext_remove(original_file)
  file_copy(original_file, copied_file, overwrite = TRUE)
  file_chmod(copied_file, "u+w")
}

# common functions and variables

MUT_FILE <- "pnet_prostate_paper/_database/prostate/processed/P1000_final_analysis_set_cross_important_only.csv"
CNV_FILE <- "pnet_prostate_paper/_database/prostate/processed/P1000_data_CNA_paper.csv"
LABEL_FILE <- "pnet_prostate_paper/_database/prostate/processed/response_paper.csv"

MUT_FILE_ORIGINAL <- "pnet_data/P1000_final_analysis_set_cross_important_only.csv"
CNV_FILE_ORIGINAL <- "pnet_data/P1000_data_CNA_paper.csv"
LABEL_FILE_ORIGINAL <- "pnet_data/response_paper.csv"

restore_input_files <- function() {
  walk2(
    c(MUT_FILE_ORIGINAL, CNV_FILE_ORIGINAL, LABEL_FILE_ORIGINAL),
    c(MUT_FILE, CNV_FILE, LABEL_FILE),
    file_copy,
    overwrite = TRUE
  )
}

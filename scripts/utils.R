
# Constants ---------------------------------------------------------------

# files that will be mounted into the PNET Docker container
MOUNTED_FILES <- list(
  mutations = "pnet_data/mounted/mutations.csv",
  cnvs = "pnet_data/mounted/cnvs.csv",
  labels = "pnet_data/mounted/labels.csv",
  training_set = "pnet_data/mounted/training_set.csv",
  test_set = "pnet_data/mounted/test_set.csv",
  validation_set = "pnet_data/mounted/validation_set.csv"
)

# original PNET input files from _database.zip
ORIGINAL_FILES <- list(
  mutations = "pnet_data/original/P1000_final_analysis_set_cross_important_only.csv",
  cnvs = "pnet_data/original/P1000_data_CNA_paper.csv",
  labels = "pnet_data/original/response_paper.csv",
  training_set = "pnet_data/original/training_set_0.csv",
  test_set = "pnet_data/original/test_set.csv",
  validation_set = "pnet_data/original/validation_set.csv"
)



# Logging -----------------------------------------------------------------

default_logger <- log4r::logger(threshold = "DEBUG")

debug <- function(..., .envir = parent.frame()) {
  log4r::debug(default_logger, glue::glue(..., .envir = .envir))
}

info <- function(..., .envir = parent.frame()) {
  log4r::info(default_logger, glue::glue(..., .envir = .envir))
}

warn <- function(..., .envir = parent.frame()) {
  log4r::warn(default_logger, glue::glue(..., .envir = .envir))
}


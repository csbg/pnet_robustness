# load original PNET data files

library(tidyverse)
library(fs)
source("scripts/utils.R")


dir_create("pnet_data/mounted")
walk2(
  ORIGINAL_FILES,
  MOUNTED_FILES,
  file_copy,
  overwrite = TRUE
)

info("Loaded original PNET input data")

library(tidyverse)
library(fs)
source("common_functions.R")


restore_input_files()


samples_metastatic <-
  read_csv(LABEL_FILE_ORIGINAL, col_types = "ci") %>%
  filter(response == 1L) %>%
  pull(id)


mutations <- read_csv(MUT_FILE_ORIGINAL)

mutations %>%
  mutate(
    across(!1, ~as.integer(Tumor_Sample_Barcode %in% samples_metastatic))
  ) %>%
  write_csv(MUT_FILE)



cnvs <- read_csv(CNV_FILE_ORIGINAL)

c("", colnames(cnvs)[-1]) %>%
  str_c(collapse = ",") %>%
  write_lines(CNV_FILE)

cnvs %>%
  mutate(
    across(!1, ~as.integer(X1 %in% samples_metastatic) * 2)
  ) %>%
  write_csv(CNV_FILE, append = TRUE)

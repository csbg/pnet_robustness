library(tidyverse)
library(fs)
source("common_functions.R")


restore_input_file(LABEL_FILE)


samples_metastatic <-
  read_csv(LABEL_FILE, col_types = "ci") %>%
  filter(response == 1L) %>%
  pull(id)


mutations <- read_csv(MUT_FILE)

mutations %>%
  mutate(
    across(!1, ~as.integer(Tumor_Sample_Barcode %in% samples_metastatic))
  ) %>%
  write_csv(path_ext_remove(MUT_FILE))



cnvs <- read_csv(CNV_FILE)

c("", colnames(cnvs)[-1]) %>%
  str_c(collapse = ",") %>%
  write_lines(path_ext_remove(CNV_FILE))

cnvs %>%
  mutate(
    across(!1, ~as.integer(X1 %in% samples_metastatic) * 2)
  ) %>%
  write_csv(path_ext_remove(CNV_FILE), append = TRUE)

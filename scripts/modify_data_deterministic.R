# make mutation and CNV data deterministic

library(tidyverse)
library(fs)
source("scripts/utils.R")


samples_metastatic <-
  read_csv(MOUNTED_FILES$labels, col_types = "ci") %>%
  filter(response == 1L) %>%
  pull(id)


mutations <- read_csv(MOUNTED_FILES$mutations)

mutations %>%
  mutate(
    across(!1, ~as.integer(Tumor_Sample_Barcode %in% samples_metastatic))
  ) %>%
  write_csv(MOUNTED_FILES$mutations)



cnvs <-
  read_csv(MOUNTED_FILES$cnvs) %>%
  mutate(across(!1, ~as.integer(...1 %in% samples_metastatic) * 2))

colnames(cnvs)[1] <- ""

cnvs %>% write_csv(MOUNTED_FILES$cnvs)

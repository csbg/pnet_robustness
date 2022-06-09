# scramble input feature
# may process two command line arguments:
# (1) probability for a nonzero feature value
# (2) random seed

library(tidyverse)
library(fs)
source("common_functions.R")



# Parameters --------------------------------------------------------------

# process command line arguments
# if none is supplied or interactive mode, keep the class frequency
args <- commandArgs(TRUE)
if (length(args) == 0L) {
  PROB_MUT_CNV <- 0.5
  SEED <- 0
} else {
  PROB_MUT_CNV <- as.numeric(args[1])
  SEED <- as.numeric(args[2])
}



# Change ------------------------------------------------------------------

set.seed(SEED)

# make 1000 samples with equally distributed labels
labels <- read_csv(LABEL_FILE)

labels %>%
  mutate(response = rep_len(0:1, nrow(labels))) %>%
  write_csv(path_ext_remove(LABEL_FILE))


# assign random mutations (original gene names)
mut_colnames <-
  read_csv(MUT_FILE) %>%
  colnames()

bind_cols(
  select(labels, !!mut_colnames[1] := id),
  map_dfc(
    mut_colnames[-1],
    ~tibble(
      !!.x := sample(
        c(0L, 1L),
        nrow(labels),
        replace = TRUE,
        prob = c(1 - PROB_MUT_CNV, PROB_MUT_CNV)
      )
    )
  )
) %>%
  write_csv(path_ext_remove(MUT_FILE))


# assign random CNVs (original gene names)
cnv_colnames <-
  read_csv(CNV_FILE) %>%
  colnames()

c("", cnv_colnames[-1]) %>%
  str_c(collapse = ",") %>%
  write_lines(path_ext_remove(CNV_FILE))

bind_cols(
  select(labels, id),
  map_dfc(
    cnv_colnames[-1],
    ~tibble(
      !!.x := sample(
        c(0L, 2L),
        nrow(labels),
        replace = TRUE,
        prob = c(1 - PROB_MUT_CNV, PROB_MUT_CNV)
      )
    )
  )
) %>%
  write_csv(path_ext_remove(CNV_FILE), append = TRUE)

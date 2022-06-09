# scramble input feature
# may process two command line arguments:
# (1) probability for a nonzero feature value
# (2) random seed

library(tidyverse)
library(fs)
source("common_functions.R")


restore_input_files()

# process command line arguments
# if none is supplied or interactive mode, keep the class frequency
args <- commandArgs(TRUE)
if (length(args) == 0L) {
  prob_mut_cnv <- 0.5
  seed <- 0
} else {
  prob_mut_cnv <- as.numeric(args[1])
  seed <- as.numeric(args[2])
}


set.seed(seed)

# make 1000 samples with equally distributed labels
labels <- read_csv(LABEL_FILE_ORIGINAL)

labels %>%
  mutate(response = rep_len(0:1, nrow(labels))) %>%
  write_csv(LABEL_FILE)


# assign random mutations (original gene names)
mut_colnames <-
  read_csv(MUT_FILE_ORIGINAL) %>%
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
        prob = c(1 - prob_mut_cnv, prob_mut_cnv)
      )
    )
  )
) %>%
  write_csv(MUT_FILE)


# assign random CNVs (original gene names)
cnv_colnames <-
  read_csv(CNV_FILE_ORIGINAL) %>%
  colnames()

c("", cnv_colnames[-1]) %>%
  str_c(collapse = ",") %>%
  write_lines(CNV_FILE)

bind_cols(
  select(labels, id),
  map_dfc(
    cnv_colnames[-1],
    ~tibble(
      !!.x := sample(
        c(0L, 2L),
        nrow(labels),
        replace = TRUE,
        prob = c(1 - prob_mut_cnv, prob_mut_cnv)
      )
    )
  )
) %>%
  write_csv(CNV_FILE, append = TRUE)

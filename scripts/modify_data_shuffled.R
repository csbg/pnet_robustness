# shuffle input labels
#
# the script has two optional command line arguments:
# (1) a Boolean value that decides whether to preserve the class frequency
# (2) random seed for scrambling

library(tidyverse)
library(fs)
source("scripts/utils.R")


# process command line arguments
args <- commandArgs(TRUE)
if (length(args) == 0L) {
  keep_class_frequency <- FALSE
  seed <- 0
} else {
  keep_class_frequency <- as.logical(args[1])
  seed <- as.numeric(args[2])
}

# shuffle labels
labels <- read_csv(MOUNTED_FILES$labels)
label_1_frq <- if_else(keep_class_frequency, mean(labels$response), 0.5)

set.seed(seed)
shuffled_response <- sample(
  0:1,
  nrow(labels),
  replace = TRUE,
  prob = c(1 - label_1_frq, label_1_frq)
)

labels %>%
  mutate(response = shuffled_response) %>%
  write_csv(MOUNTED_FILES$labels)

info("Shuffled {nrow(labels)} labels")

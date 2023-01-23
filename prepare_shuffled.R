# scramble input labels in response_paper.csv
# optionally ensure that class frequency remains constant
# can process two command line argument:
# (1) a Boolean value that decides whether to preserve the class frequency
# (2) random seed for scrambling

library(tidyverse)
library(fs)
source("common_functions.R")


restore_input_files()

# process command line arguments
args <- commandArgs(TRUE)
if (length(args) == 0L) {
  keep_class_frequency <- FALSE
  seed <- 0
} else {
  keep_class_frequency <- as.logical(args[1])
  seed <- as.numeric(args[2])
}

# scramble labels
labels <- read_csv(LABEL_FILE_ORIGINAL)

if (keep_class_frequency) {
  label_1_frq <- mean(labels$response)
} else {
  label_1_frq <- 0.5
}

set.seed(seed)
scrambled_response <- sample(
  0:1,
  nrow(labels),
  replace = TRUE,
  prob = c(1 - label_1_frq, label_1_frq)
)

labels %>%
  mutate(response = scrambled_response) %>%
  write_csv(LABEL_FILE)

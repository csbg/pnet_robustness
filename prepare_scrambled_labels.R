# scramble input labels in response_paper.csv
# optionally ensure that class frequency remains constant
# can process two command line argument:
# (1) a Boolean value that decides whether to preserve the class frequency
# (2) random seed for scrambling

library(tidyverse)
library(fs)
source("common_functions.R")



# process command line arguments
args <- commandArgs(TRUE)
if (length(args) == 0L) {
  KEEP_CLASS_FREQUENCY <- TRUE
  SEED <- 0
} else {
  KEEP_CLASS_FREQUENCY <- as.logical(args[1])
  SEED <- as.numeric(args[2])
}

# restore unused input files
restore_input_file(MUT_FILE)
restore_input_file(CNV_FILE)


# scramble labels
scramble_labels <- function(keep_class_frequency = TRUE) {
  labels <- read_csv(LABEL_FILE)

  if (keep_class_frequency)
    label_1_frq <- mean(labels$response)
  else
    label_1_frq <- 0.5

  set.seed(SEED)
  scrambled_response <- sample(
    0:1,
    nrow(labels),
    replace = TRUE,
    prob = c(1 - label_1_frq, label_1_frq)
  )

  labels %>%
    mutate(response = scrambled_response) %>%
    write_csv(path_ext_remove(LABEL_FILE))
}

scramble_labels(KEEP_CLASS_FREQUENCY)

# simply copies the .csv.original files to .csv

library(tidyverse)
library(fs)
source("common_functions.R")

restore_input_file(LABEL_FILE)
restore_input_file(MUT_FILE)
restore_input_file(CNV_FILE)

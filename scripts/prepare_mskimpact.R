# load data from the MSK-IMPACT dataset

library(tidyverse)
library(fs)
source("scripts/common_functions.R")



# Load data ---------------------------------------------------------------

# labels for 10,945 samples
label_data <-
  read_tsv("pnet_data/msk_impact_2017/data_clinical_sample.txt", skip = 4) %>%
  mutate(
    .keep = "none",
    id = SAMPLE_ID,
    response = if_else(SAMPLE_TYPE == "Primary", 0L, 1L),
  )

# mutation data
mut_data <-
  read_tsv("pnet_data/msk_impact_2017/data_mutations.txt", skip = 1) %>%
  summarise(
    .by = c(Tumor_Sample_Barcode, Hugo_Symbol),
    mut = n()
  ) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = mut)

# cnv data
Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)
cnv_data <-
  read_tsv("pnet_data/msk_impact_2017/data_cna.txt") %>%
  pivot_longer(!Hugo_Symbol, names_to = "sample", values_to = "cna") %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = cna)

colnames(cnv_data)[1] <- ""

# splits
set.seed(0)
label_rows <- sample(nrow(label_data))
train_size <- as.integer(length(label_rows) * .8)
test_size <- as.integer(length(label_rows) * .1)

train_rows <- label_rows[1:train_size]
test_rows <- label_rows[(train_size + 1):(train_size + test_size)]
valid_rows <- label_rows[(train_size + test_size + 1):nrow(label_data)]

train_data <- label_data[train_rows, ]
test_data <- label_data[test_rows, ]
valid_data <- label_data[valid_rows, ]



# Save data ---------------------------------------------------------------

label_data %>% write_csv(MOUNTED_FILES$labels)
mut_data %>% write_csv(MOUNTED_FILES$mutations)
cnv_data %>% write_csv(MOUNTED_FILES$cnvs)


write_csv_with_row_index <- function(df, file) {
  df <-
    df %>%
    mutate(rn = row_number() - 1L, .before = 1)
  colnames(df)[1] <- ""
  write_csv(df, file)
}

train_data %>% write_csv_with_row_index(MOUNTED_FILES$training_set)
test_data %>% write_csv_with_row_index(MOUNTED_FILES$test_set)
valid_data %>% write_csv_with_row_index(MOUNTED_FILES$validation_set)

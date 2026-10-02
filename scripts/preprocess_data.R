# Preprocess GTEx data for the HOHC Shiny app.
#
# The app only ever uses the per-subject GTEx table through
#   group_by(Description, SMTS) %>% summarise(TPM = mean(TPM))   (within one SEX)
# so we precompute that mean once here. After aggregation every
# (Description, SMTS, SEX) group has exactly one row, which makes the
# original summarise() in sankey.R an identity operation — results are
# mathematically unchanged while the data shrinks from 156 MB to < 1 MB.
#
# Run from the project root:
#   Rscript scripts/preprocess_data.R

library(readr)
library(dplyr)

message("Reading data/GTEx_HR_combined.tsv ...")
gtex_raw <- read_tsv("data/GTEx_HR_combined.tsv", show_col_types = FALSE)

message("Aggregating mean TPM per (Description, SMTS, SEX) ...")
gtex_mean <- gtex_raw %>%
  group_by(Description, SMTS, SEX) %>%
  summarise(TPM = mean(TPM), .groups = "drop")

message("Rows before: ", nrow(gtex_raw), "  after: ", nrow(gtex_mean))

out <- "data/GTEx_HR_mean.rds"
saveRDS(gtex_mean, out, compress = "xz")
message("Saved ", out, " (", round(file.size(out) / 1024), " KB)")

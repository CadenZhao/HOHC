# Verify that the pre-aggregated data gives results identical to the
# original per-subject data for the app's actual computations.
#
# Checks, for both sexes:
#   1. The summarised table produced inside sankey.R from raw data equals
#      the one produced from the pre-aggregated data.
#   2. filter_data() (the app's download output) is identical under the
#      default UI settings and a few other filter combinations.
#
# Run from the project root AFTER scripts/preprocess_data.R:
#   Rscript scripts/verify_equivalence.R

suppressPackageStartupMessages(library(tidyverse))

gtex_raw  <- read_tsv("data/GTEx_HR_combined.tsv", show_col_types = FALSE)
gtex_mean <- readRDS("data/GTEx_HR_mean.rds")

df <- readxl::read_excel("data/hormone_receptor_list.xlsx")
df$Hormone_name1 <- recode(df$Hormone_name1, "orphan" = "Unknown")
df$Hormone_chemical_classes <- recode(df$Hormone_chemical_classes, "NA" = "Unknown")
df$Hormone1_tissue <- recode(df$Hormone1_tissue, "NA" = "Unknown")
df <- df %>% separate_rows(Hormone1_tissue, sep = ", ")
df <- df %>% filter(Hormone1_tissue != "Unknown")
vec_receptor <- df$Gene_symbol

options(dplyr.summarise.inform = FALSE)

summarise_as_app <- function(gtex, sex) {
  gtex %>%
    filter(Description %in% vec_receptor, SEX == sex) %>%
    group_by(Description, SMTS) %>%
    summarise(TPM = mean(TPM)) %>%
    ungroup() %>%
    arrange(Description, SMTS)
}

ok <- TRUE
for (sex in c("Male", "Female")) {
  a <- summarise_as_app(gtex_raw, sex)
  b <- summarise_as_app(gtex_mean, sex)
  same <- isTRUE(all.equal(as.data.frame(a), as.data.frame(b)))
  cat(sprintf("[1] summarised table identical (%s): %s\n", sex, same))
  if (!same) ok <- FALSE
}

# vec_tissue used by the UI must be unchanged
same_tissue <- identical(sort(unique(gtex_raw$SMTS)), sort(unique(gtex_mean$SMTS)))
cat(sprintf("[2] tissue vector identical: %s\n", same_tissue))
if (!same_tissue) ok <- FALSE

# Full filter_data() comparison across several filter settings
source("sankey.R")
settings <- list(
  list(sex = "Male",   hc = "All",     hn = c("Androgen", "Glucagon-Like Peptide 1"), rc = "All",        n = 3,  lo = 0.1, hi = 6.35),
  list(sex = "Female", hc = "All",     hn = c("Androgen", "Glucagon-Like Peptide 1"), rc = "All",        n = 3,  lo = 0.1, hi = 6.35),
  list(sex = "Male",   hc = "Peptide", hn = NULL,                                     rc = "Nuclear HR", n = 5,  lo = 0.5, hi = 6.0),
  list(sex = "Female", hc = "All",     hn = NULL,                                     rc = "All",        n = 10, lo = 0.0, hi = 6.35)
)
for (i in seq_along(settings)) {
  s <- settings[[i]]
  run <- function(gtex) {
    filter_data(s$sex, s$hc, s$hn, NULL, s$rc, NULL, NULL, s$n, s$lo, s$hi,
                gtex, df, vec_receptor) %>%
      arrange(across(everything()))
  }
  a <- run(gtex_raw)
  b <- run(gtex_mean)
  same <- isTRUE(all.equal(as.data.frame(a), as.data.frame(b)))
  cat(sprintf("[3.%d] filter_data identical (%s, n=%d): %s (%d rows)\n",
              i, s$sex, s$n, same, nrow(a)))
  if (!same) ok <- FALSE
}

cat(if (ok) "\nALL CHECKS PASSED\n" else "\nMISMATCH FOUND\n")
if (!ok) quit(status = 1)

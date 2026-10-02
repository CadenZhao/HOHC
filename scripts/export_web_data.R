# Export the data needed by the static web version of HOHC.
#
# Produces site/data/hohc_data.json containing:
#   - receptors: the processed hormone-receptor annotation table
#     (same processing as in app.R)
#   - expression: pre-aggregated mean TPM per (gene, organ, sex)
#   - tissues: organ vector for the UI selectors
#
# Run from the project root:
#   Rscript scripts/export_web_data.R

suppressPackageStartupMessages(library(tidyverse))

gtex_mean <- readRDS("data/GTEx_HR_mean.rds")

# Same annotation processing as app.R
df <- readxl::read_excel("data/hormone_receptor_list.xlsx")
df$Hormone_name1 <- recode(df$Hormone_name1, "orphan" = "Unknown")
df$Hormone_chemical_classes <- recode(df$Hormone_chemical_classes, "NA" = "Unknown")
df$Hormone1_tissue <- recode(df$Hormone1_tissue, "NA" = "Unknown")
df <- df %>% separate_rows(Hormone1_tissue, sep = ", ")
df <- df %>% filter(Hormone1_tissue != "Unknown")

receptors <- df %>%
  transmute(
    gene = Gene_symbol,
    hormone = Hormone_name1,
    hormone_class = Hormone_chemical_classes,
    hormone_organ = Hormone1_tissue,
    receptor_class = Receptor_subcellular
  )

expression <- gtex_mean %>%
  transmute(
    gene = Description,
    organ = SMTS,
    sex = SEX,
    tpm = round(TPM, 6)
  )

out <- list(
  receptors = receptors,
  expression = expression,
  tissues = sort(unique(gtex_mean$SMTS))
)

dir.create("site/data", recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(out, "site/data/hohc_data.json", auto_unbox = FALSE, digits = NA)
message("Wrote site/data/hohc_data.json (",
        round(file.size("site/data/hohc_data.json") / 1024), " KB)")

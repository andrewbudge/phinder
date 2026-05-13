#!/usr/bin/env Rscript
# Filter geNomad virus_summary output.
#
# Rules:
#   * topology == "DTR"                                   -> keep regardless of score
#   * topology == "Provirus" & virus_score >= MIN_PROV    -> keep
#   * everything else                                     -> drop
#
# Usage: filter_genomad.R <input_summary.tsv> <output.tsv> <min_provirus_score>

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: filter_genomad.R <input_summary.tsv> <output.tsv> <min_provirus_score>")
}
input_tsv    <- args[1]
output_tsv   <- args[2]
min_prov_sc  <- as.numeric(args[3])

genomad <- read_tsv(input_tsv, show_col_types = FALSE)

cat("Input rows: ", nrow(genomad), "\n", sep = "")
cat("Topology breakdown:\n")
print(table(genomad$topology, useNA = "ifany"))

filtered <- genomad |>
  filter(topology == "DTR" |
         (topology == "Provirus" & virus_score >= min_prov_sc))

cat("\nKept rows: ", nrow(filtered), "\n", sep = "")
cat("Kept topology breakdown:\n")
print(table(filtered$topology, useNA = "ifany"))

write_tsv(filtered, output_tsv)
cat("\nWrote: ", output_tsv, "\n", sep = "")

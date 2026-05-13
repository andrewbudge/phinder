#!/usr/bin/env Rscript
# Filter CheckV quality_summary and join with the geNomad filter to produce
# the candidate-phage table that drives Pharokka / PhaBOX.
#
# Rules:
#   * Keep if checkv_quality %in% quality_keep  OR  topology == "DTR"
#       - DTR bypass: CheckV under-scores DTRs because they lack host flanking
#   * Flag low coverage (< min_coverage) but do NOT drop
#       - Coverage is parsed opportunistically from the input FASTA headers
#         (assemblers like metaMDBG embed `coverage=N` in headers). If the
#         header lacks coverage info, the flag column is left as NA.
#
# Usage:
#   filter_checkv.R <filtered_genomad.tsv> <checkv_summary.tsv> \
#                   <input_contigs.fa[.gz]> <output.tsv> \
#                   <min_coverage> <quality_keep_csv>

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: filter_checkv.R <filtered_genomad.tsv> <checkv_summary.tsv> ",
       "<input_contigs.fa[.gz]> <output.tsv> <min_coverage> <quality_keep_csv>")
}
genomad_path  <- args[1]
checkv_path   <- args[2]
contigs_path  <- args[3]
out_path      <- args[4]
min_coverage  <- as.numeric(args[5])
quality_keep  <- strsplit(args[6], ",", fixed = TRUE)[[1]]

filtered_genomad <- read_tsv(genomad_path, show_col_types = FALSE)
checkv           <- read_tsv(checkv_path,  show_col_types = FALSE)

# --- Coverage map (best-effort) ---------------------------------------------
read_headers <- function(f) {
  con <- if (grepl("\\.gz$", f)) gzfile(f, "rt") else file(f, "rt")
  on.exit(close(con))
  lines <- grep("^>", readLines(con), value = TRUE)
  tibble(header = sub("^>", "", lines))
}

headers <- read_headers(contigs_path)

coverage_map <- headers |>
  mutate(
    base_id  = str_extract(header, "^\\S+"),
    coverage = suppressWarnings(
      as.numeric(str_match(header, "coverage=([0-9.]+)")[, 2])
    )
  ) |>
  select(base_id, coverage)

has_coverage <- any(!is.na(coverage_map$coverage))

attach_coverage <- function(df, id_col) {
  df |>
    mutate(base_id = sub("\\|provirus_.*$", "", .data[[id_col]])) |>
    left_join(coverage_map, by = "base_id") |>
    select(-base_id)
}

# --- Join + filter ----------------------------------------------------------
potential_phages <- filtered_genomad |>
  left_join(checkv, by = c("seq_name" = "contig_id")) |>
  attach_coverage("seq_name") |>
  filter(checkv_quality %in% quality_keep | topology == "DTR") |>
  mutate(low_coverage = !is.na(coverage) & coverage < min_coverage) |>
  select(seq_name, length, topology, virus_score, n_hallmarks,
         completeness, checkv_quality, coverage, low_coverage, taxonomy)

cat("Candidate phages: ", nrow(potential_phages), "\n", sep = "")
if (has_coverage) {
  cat("Flagged low coverage (<", min_coverage, "x): ",
      sum(potential_phages$low_coverage, na.rm = TRUE), "\n", sep = "")
} else {
  cat("Note: no `coverage=` field found in input FASTA headers; ",
      "low_coverage flagging skipped.\n", sep = "")
}
cat("Quality breakdown:\n")
print(table(potential_phages$checkv_quality, useNA = "ifany"))

write_tsv(potential_phages, out_path)
cat("\nWrote: ", out_path, "\n", sep = "")

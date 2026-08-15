# Phase 6: exports - cleaned dataset, QC flags, and a text/markdown QC report.
# Pure functions (no Shiny).

library(tidyverse)

export_clean_long <- function(data, path) {
  readr::write_csv(data, path)
}

export_qc_flags <- function(flags, path) {
  readr::write_csv(flags, path)
}

#' Build a Markdown QC report string: equipment metadata (from the raw CSV),
#' MFI thresholds, bead-count QC summary, and (if supplied) the removal log
#' from apply_qc_removals().
build_qc_report <- function(parsed, thresholds, flags, removal_log = NULL) {
  meta <- parsed$metadata
  meta_lines <- sprintf("- **%s**: %s", meta$field, ifelse(is.na(meta$value), "", meta$value))

  bead_summary <- flags |>
    dplyr::count(bead_quality) |>
    (\(x) sprintf("- %s: %d wells", x$bead_quality, x$n))()

  review_n <- nrow(flagged_for_review(flags))

  lines <- c(
    "# Serology QC Report",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Equipment / batch metadata",
    meta_lines,
    "",
    "## MFI thresholds",
    sprintf("- Negative/background threshold (mean %s Median MFI): %.1f",
            thresholds$empty_antigen, thresholds$min),
    sprintf("- Positive/ceiling threshold (mean %s Median MFI, samples+pools): %.1f",
            thresholds$pos_antigen, thresholds$max),
    "",
    "## Bead-count QC summary (384 wells)",
    bead_summary,
    "",
    sprintf("**%d of %d wells** flagged for review (Bad/Moderate bead count or low MFI signal).",
            review_n, nrow(flags))
  )

  if (!is.null(removal_log)) {
    lines <- c(lines, "", "## Removal decisions applied", paste0("- ", removal_log))
  }

  paste(lines, collapse = "\n")
}

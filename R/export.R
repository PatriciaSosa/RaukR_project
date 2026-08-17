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

#' Export the QC report as a PDF: a plain-text summary page (same content as
#' build_qc_report(), with Markdown syntax stripped) followed by the 384-well
#' bead-count plate ("min" view) and the top 3 MFI plots as static figures.
#' Uses only base R graphics devices, no external PDF/LaTeX tools required.
export_qc_report_pdf <- function(path, parsed, thresholds, flags, merged,
                                 removal_log = NULL) {
  txt <- build_qc_report(parsed, thresholds, flags, removal_log)
  lines <- strsplit(txt, "\n")[[1]]
  lines <- gsub("\\*\\*", "", lines)
  lines <- sub("^## ", "", lines)
  lines <- sub("^# ", "", lines)

  grDevices::pdf(path, width = 8.5, height = 11)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Text summary, paginated ~50 lines per page.
  page_size <- 50
  pages <- split(lines, ceiling(seq_along(lines) / page_size))
  for (pg in pages) {
    grid::grid.newpage()
    grid::grid.text(
      paste(pg, collapse = "\n"), x = 0.04, y = 0.97, just = c("left", "top"),
      gp = grid::gpar(fontfamily = "mono", fontsize = 9)
    )
  }

  print(plot_plate_beadcount(merged, view = "min", interactive = FALSE))
  print(plot_mfi_vs_controls(merged, interactive = FALSE))
  print(plot_mfi_by_antigen(merged, interactive = FALSE))
  print(plot_mfi_by_sample(merged, thresholds = thresholds, interactive = FALSE))

  invisible(path)
}

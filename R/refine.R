# Flagging and removal of poor-quality wells/samples.
# Pure functions (no Shiny); removal happens at the well/sample level, since
# bead count and every antigen reading in a well belong to the same sample.

library(tidyverse)

# ---------------------------------------------------------------------------
# Combined per-well/sample QC flags
# ---------------------------------------------------------------------------

#' One row per well summarising bead-count quality (Phase 3). MFI-based
#' removal is always the analyst's explicit decision (made by reviewing the
#' MFI plots in the "MFI QC" tab), so no automatic MFI threshold flag is
#' computed here.
build_sample_flags <- function(merged, bead_view = "min",
                               bad_max = QC_BAD_MAX, good_min = QC_GOOD_MIN) {
  bead <- well_bead_summary(merged, view = bead_view, bad_max = bad_max, good_min = good_min) |>
    dplyr::rename(bead_quality = quality)

  bead |>
    dplyr::select(well_384, sample_id, sample_type, Quadrant,
                  bead_count, limiting_antigen = limiting_analyte, bead_quality) |>
    dplyr::mutate(
      flag_reason = dplyr::case_when(
        bead_quality == "Bad" ~ "Bad bead count",
        bead_quality == "Moderate" ~ "Moderate bead count",
        TRUE ~ "OK"
      )
    )
}

#' Subset of flags worth showing the analyst for manual review: anything
#' short of Good bead count. Mirrors the "Download QC flags" export.
flagged_for_review <- function(flags) {
  dplyr::filter(flags, bead_quality != "Good")
}

# ---------------------------------------------------------------------------
# Sample removal template (analyst decides removals manually from MFI review)
# ---------------------------------------------------------------------------

#' Every current sample_id, for the analyst to download, edit down to just
#' the sample_ids they want removed (e.g. after reviewing the MFI plots),
#' and re-upload via the manual removal list. Removal is never automatic
#' from the MFI thresholds - it is always the analyst's explicit decision.
sample_removal_template <- function(merged) {
  merged |>
    dplyr::distinct(sample_id, sample_type, well_384, Quadrant) |>
    dplyr::arrange(sample_type, sample_id)
}

# ---------------------------------------------------------------------------
# Manual removal list (analyst-identified lab errors)
# ---------------------------------------------------------------------------

#' Read an uploaded manual-removal list. Accepts a "sample_id" column
#' (case-insensitive); returns a character vector of sample_ids.
read_manual_removal_list <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  col <- names(df)[tolower(names(df)) == "sample_id"]
  if (length(col) == 0) {
    stop("Uploaded removal list must have a 'sample_id' column.")
  }
  unique(stats::na.omit(df[[col[1]]]))
}

# ---------------------------------------------------------------------------
# Apply removal decisions
# ---------------------------------------------------------------------------

#' Apply the analyst's QC decisions to the merged long table.
#'
#' remove_bead_categories: bead-count quality categories to drop, e.g.
#'   character(0) (keep all), "Bad", or c("Bad", "Moderate").
#' manual_sample_ids: additional sample_ids to drop regardless of flags
#'   (e.g. a known lab/pipetting error, or a low-MFI sample identified by
#'   the analyst while reviewing the MFI plots).
#'
#' Returns list(clean, removed, flags, log): the filtered long table, the
#' rows removed with their reasons, the full flag table, and a text log
#' suitable for the QC report.
apply_qc_removals <- function(merged, flags,
                              remove_bead_categories = character(0),
                              manual_sample_ids = character(0)) {
  decisions <- flags |>
    dplyr::mutate(
      drop_bead_count = bead_quality %in% remove_bead_categories,
      drop_manual = sample_id %in% manual_sample_ids,
      drop = drop_bead_count | drop_manual,
      removal_reason = dplyr::case_when(
        drop_manual & drop_bead_count ~ paste0("Manual removal + ", flag_reason),
        drop_manual ~ "Manual removal",
        drop_bead_count ~ flag_reason,
        TRUE ~ NA_character_
      )
    )

  dropped <- dplyr::filter(decisions, drop)
  clean <- dplyr::filter(merged, !well_384 %in% dropped$well_384)

  n_wells_total <- dplyr::n_distinct(merged$well_384)
  n_wells_dropped <- dplyr::n_distinct(dropped$well_384)

  log_lines <- c(
    sprintf("QC refinement run: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Bead-count categories removed: %s",
            if (length(remove_bead_categories)) paste(remove_bead_categories, collapse = ", ") else "none"),
    sprintf("Manual sample removals requested: %d (%s)",
            length(manual_sample_ids),
            if (length(manual_sample_ids)) paste(manual_sample_ids, collapse = ", ") else "none"),
    sprintf("Wells removed: %d of %d (%.1f%%)",
            n_wells_dropped, n_wells_total, 100 * n_wells_dropped / n_wells_total)
  )

  list(clean = clean, removed = dropped, flags = decisions, log = log_lines)
}

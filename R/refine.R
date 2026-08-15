# Flagging and removal of poor-quality wells/samples.
# Pure functions (no Shiny); removal happens at the well/sample level, since
# bead count and every antigen reading in a well belong to the same sample.

library(tidyverse)

# ---------------------------------------------------------------------------
# Combined per-well/sample QC flags
# ---------------------------------------------------------------------------

#' One row per well combining bead-count quality (Phase 3) with an MFI
#' "low signal" heuristic: the fraction of antigens whose Median MFI falls
#' below the negative/background threshold (Phase 4's `mfi_control_thresholds`
#' min). A well is flagged `low_signal_flag` when that fraction is at or above
#' `mfi_frac_threshold` (default 0.8: essentially no detectable response to
#' any antigen, suggesting a handling error rather than true seronegativity).
#'
#' This heuristic is a starting point for manual review, not an automatic
#' verdict - the analyst inspects flags before removing anything.
build_sample_flags <- function(merged, bead_view = "min",
                               bad_max = QC_BAD_MAX, good_min = QC_GOOD_MIN,
                               thresholds = NULL, mfi_frac_threshold = 0.8,
                               exclude_controls = TRUE, controls = CONTROL_LABELS) {
  bead <- well_bead_summary(merged, view = bead_view, bad_max = bad_max, good_min = good_min) |>
    dplyr::rename(bead_quality = quality)

  mfi_df <- merged
  if (exclude_controls) mfi_df <- dplyr::filter(mfi_df, !antigen %in% controls)
  if (is.null(thresholds)) thresholds <- mfi_control_thresholds(merged)

  mfi_flags <- mfi_df |>
    dplyr::group_by(well_384) |>
    dplyr::summarise(
      n_antigens = dplyr::n(),
      n_below_min = sum(median_mfi < thresholds$min, na.rm = TRUE),
      frac_below_min = n_below_min / n_antigens,
      .groups = "drop"
    ) |>
    dplyr::mutate(low_signal_flag = frac_below_min >= mfi_frac_threshold)

  bead |>
    dplyr::select(well_384, sample_id, sample_type, Quadrant,
                  bead_count, limiting_antigen = limiting_analyte, bead_quality) |>
    dplyr::left_join(mfi_flags, by = "well_384") |>
    dplyr::mutate(
      flag_reason = dplyr::case_when(
        bead_quality == "Bad" & low_signal_flag ~ "Bad bead count; low MFI signal",
        bead_quality == "Bad" ~ "Bad bead count",
        bead_quality == "Moderate" & low_signal_flag ~ "Moderate bead count; low MFI signal",
        bead_quality == "Moderate" ~ "Moderate bead count",
        low_signal_flag ~ "Low MFI signal",
        TRUE ~ "OK"
      )
    )
}

#' Subset of flags worth showing the analyst for manual review: anything not
#' both Good bead count and MFI-clean. Mirrors the "Download QC flags" export.
flagged_for_review <- function(flags) {
  dplyr::filter(flags, bead_quality != "Good" | low_signal_flag)
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
#' remove_low_mfi: drop wells flagged `low_signal_flag`.
#' manual_sample_ids: additional sample_ids to drop regardless of flags
#'   (e.g. a known lab/pipetting error).
#'
#' Returns list(clean, removed, flags, log): the filtered long table, the
#' rows removed with their reasons, the full flag table, and a text log
#' suitable for the QC report.
apply_qc_removals <- function(merged, flags,
                              remove_bead_categories = character(0),
                              remove_low_mfi = FALSE,
                              manual_sample_ids = character(0)) {
  decisions <- flags |>
    dplyr::mutate(
      drop_bead_count = bead_quality %in% remove_bead_categories,
      drop_low_mfi = remove_low_mfi & low_signal_flag,
      drop_manual = sample_id %in% manual_sample_ids,
      drop = drop_bead_count | drop_low_mfi | drop_manual,
      removal_reason = dplyr::case_when(
        drop_manual & (drop_bead_count | drop_low_mfi) ~
          paste0("Manual removal + ", flag_reason),
        drop_manual ~ "Manual removal",
        drop_bead_count | drop_low_mfi ~ flag_reason,
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
    sprintf("Low-MFI-signal removal applied: %s", remove_low_mfi),
    sprintf("Manual sample removals requested: %d (%s)",
            length(manual_sample_ids),
            if (length(manual_sample_ids)) paste(manual_sample_ids, collapse = ", ") else "none"),
    sprintf("Wells removed: %d of %d (%.1f%%)",
            n_wells_dropped, n_wells_total, 100 * n_wells_dropped / n_wells_total)
  )

  list(clean = clean, removed = dropped, flags = decisions, log = log_lines)
}

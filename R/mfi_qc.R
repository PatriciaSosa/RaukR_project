# Median MFI QC: control-anchored thresholds, abundance plots, and pool CV%.
# Pure functions (no Shiny); interactive plots are plotly widgets.

library(tidyverse)

# Default control-bead names, matching the assay design in the project plan.
# Override these if a real antigen table uses different labels.
CONTROL_LABELS <- c(empty = "Bare", pos_igg = "ahIgG",
                    neg_igm = "ahIgM", semi_pos = "EBNA1")

is_control_antigen <- function(antigen, controls = CONTROL_LABELS) {
  antigen %in% controls
}

#' Identify control-bead antigen names from the antigen table's own metadata
#' instead of hardcoding literal antigen names. Antigen names (and even the
#' virus panel itself) differ from study to study, but the *role* of the
#' control beads is recorded consistently in Antigen_Group:
#'   - "Empty_bead"        -> the empty/bare bead
#'   - "Anti-human_IgG"    -> positive control bead
#'   - "Anti-human_IgM"    -> negative control bead
#' The semi-positive control (e.g. EBNA1) is an analyst's assay-design
#' choice, not something flagged in Antigen_Group/virus_genus, so it can't
#' be derived this way; it falls back to `defaults["semi_pos"]` (or NA if
#' that antigen isn't present in this particular panel).
derive_control_labels <- function(antigen_tbl, defaults = CONTROL_LABELS) {
  pick <- function(group_name) {
    hit <- antigen_tbl$antigen[trimws(antigen_tbl$Antigen_Group) == group_name]
    if (length(hit) == 0 || is.na(hit[1])) NA_character_ else hit[1]
  }
  labels <- c(
    empty    = pick("Empty_bead"),
    pos_igg  = pick("Anti-human_IgG"),
    neg_igm  = pick("Anti-human_IgM"),
    semi_pos = NA_character_
  )
  # Fall back to the built-in defaults for anything not found via metadata
  # (older antigen tables without Antigen_Group, or semi_pos).
  needs_default <- is.na(labels) | !labels %in% antigen_tbl$antigen
  fallback <- defaults[names(labels)]
  labels[needs_default] <- fallback[needs_default]
  labels[!is.na(labels)]
}

# ---------------------------------------------------------------------------
# Median MFI of the control beads themselves (empty, ahIgG, ahIgM)
# ---------------------------------------------------------------------------

#' Scatter (jittered strip plot) of Median MFI for the control beads only
#' (empty/Bare, ahIgG, ahIgM by default): x = bead, y = MFI, colored by
#' sample type. Lets the analyst confirm each control behaves as expected
#' (ahIgG high in samples/pools, empty/ahIgM low everywhere) and spot outlier
#' wells by hovering.
plot_mfi_vs_controls <- function(merged,
                                 controls = CONTROL_LABELS[c("empty", "neg_igm", "pos_igg")],
                                 interactive = TRUE) {
  df <- merged |>
    dplyr::filter(antigen %in% controls) |>
    dplyr::mutate(
      antigen = factor(antigen, levels = controls),
      tooltip = sprintf("Sample: %s (%s)\nBead: %s\nWell: %s\nMFI: %.0f",
                        sample_id, sample_type, antigen, well_384, median_mfi)
    )

  p <- suppressWarnings(
    ggplot(df, aes(x = antigen, y = median_mfi, color = sample_type)) +
      geom_jitter(aes(text = tooltip), width = 0.15, alpha = 0.6, size = 1.4) +
      scale_color_manual(values = SAMPLE_TYPE_COLORS, name = "Sample type") +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Antigens", y = "Raw MFI Values",
           title = "Raw MFI value vs control beads") +
      theme_minimal(base_size = 11)
  )

  if (interactive) to_plotly(p) else p
}

# ---------------------------------------------------------------------------
# Boxplot: Median MFI by sample (identify low/high samples by clinical code)
# ---------------------------------------------------------------------------

#' Boxplot of Median MFI per sample (one box per sample_id across its
#' antigens), colored by sample type. Hovering a point reveals the clinical
#' code. Optionally restrict to a Quadrant/sample_type subset for readability.
plot_mfi_by_sample <- function(merged, quadrant = NULL, sample_types = NULL,
                               exclude_controls = TRUE, controls = CONTROL_LABELS,
                               interactive = TRUE) {
  df <- merged
  if (exclude_controls) df <- dplyr::filter(df, !antigen %in% controls)
  if (!is.null(quadrant)) df <- dplyr::filter(df, Quadrant %in% quadrant)
  if (!is.null(sample_types)) df <- dplyr::filter(df, sample_type %in% sample_types)

  ord <- df |> dplyr::distinct(sample_id, sample_type) |>
    dplyr::arrange(sample_type, sample_id) |> dplyr::pull(sample_id)
  df <- df |>
    dplyr::mutate(
      sample_id = factor(sample_id, levels = ord),
      tooltip = sprintf("Sample: %s (%s)\nAntigen: %s\nMFI: %.0f",
                        sample_id, sample_type, antigen, median_mfi)
    )

  p <- suppressWarnings(
    ggplot(df, aes(x = sample_id, y = median_mfi)) +
      geom_boxplot(aes(color = sample_type), outlier.shape = NA, fill = NA, linewidth = 0.3) +
      geom_jitter(aes(color = sample_type, text = tooltip), width = 0.2, alpha = 0.5, size = 0.9) +
      scale_color_manual(values = SAMPLE_TYPE_COLORS, name = "Sample type") +
      scale_y_log10() +
      labs(x = "Cohort", y = "Signal (log scale)",
           title = "Raw MFI value vs all samples",
           subtitle = "Signals across samples") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  )

  if (interactive) to_plotly(p) else p
}

# ---------------------------------------------------------------------------
# Boxplot: Median MFI by antigen (samples overlaid, colored by sample type)
# ---------------------------------------------------------------------------

#' Boxplot of raw Median MFI per antigen (all analytes, including controls,
#' by default) with individual sample points overlaid. Boxes are colored by
#' Antigen_Group (pathogen/role) and points by sample_type - two separate
#' legends, matching the reference QC figures.
plot_mfi_by_antigen <- function(merged, exclude_controls = FALSE,
                                controls = CONTROL_LABELS, interactive = TRUE) {
  df <- merged
  if (exclude_controls) df <- dplyr::filter(df, !antigen %in% controls)
  ord <- df |> dplyr::distinct(antigen, Antigen_Group) |>
    dplyr::arrange(Antigen_Group, antigen) |> dplyr::pull(antigen)
  df <- df |>
    dplyr::mutate(
      antigen = factor(antigen, levels = ord),
      tooltip = sprintf("Sample: %s (%s)\nAntigen: %s\nMFI: %.0f",
                        sample_id, sample_type, antigen, median_mfi)
    )

  p <- suppressWarnings(
    ggplot(df, aes(x = antigen, y = median_mfi)) +
      geom_boxplot(aes(fill = Antigen_Group), outlier.shape = NA, alpha = 0.3, linewidth = 0.3) +
      geom_jitter(aes(color = sample_type, text = tooltip), width = 0.15, alpha = 0.5, size = 0.8) +
      scale_color_manual(values = SAMPLE_TYPE_COLORS, name = "Sample type") +
      scale_y_continuous(labels = scales::comma) +
      guides(fill = guide_legend(title = "Antigen group")) +
      labs(x = "Antigens", y = "Raw MFI Values",
           title = "Raw MFI value vs antigens") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  if (interactive) to_plotly(p) else p
}

# ---------------------------------------------------------------------------
# Boxplot: Median MFI by sample_type, faceted by antigen (user-selectable)
# ---------------------------------------------------------------------------

#' Median MFI vs sample type for one (or more) chosen antigens. x-axis order
#' is fixed to blank -> pool -> sample regardless of alphabetical/factor order.
plot_mfi_by_antigen_facet <- function(merged, antigens = NULL, interactive = TRUE) {
  df <- merged
  if (!is.null(antigens)) df <- dplyr::filter(df, antigen %in% antigens)
  df <- df |>
    dplyr::mutate(
      sample_type = factor(sample_type, levels = c("blank", "pool", "sample")),
      tooltip = sprintf("Sample: %s\nMFI: %.0f", sample_id, median_mfi)
    )

  p <- suppressWarnings(
    ggplot(df, aes(x = sample_type, y = median_mfi, color = sample_type)) +
      geom_boxplot(outlier.shape = NA, fill = NA) +
      geom_jitter(aes(text = tooltip), width = 0.2, alpha = 0.4, size = 1.2) +
      scale_color_manual(values = SAMPLE_TYPE_COLORS, guide = "none") +
      scale_y_log10() +
      facet_wrap(~antigen, scales = "free_y") +
      labs(x = NULL, y = "Median MFI (log scale)",
           title = "Median MFI by sample type") +
      theme_minimal(base_size = 12)
  )

  if (interactive) to_plotly(p) else p
}

# ---------------------------------------------------------------------------
# CV% by antigen, pools only (assay reproducibility)
# ---------------------------------------------------------------------------

compute_pool_cv <- function(merged, pool_type = "pool",
                            exclude_controls = TRUE, controls = CONTROL_LABELS) {
  df <- merged
  if (exclude_controls) df <- dplyr::filter(df, !antigen %in% controls)
  df |>
    dplyr::filter(sample_type == pool_type) |>
    dplyr::group_by(antigen) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_mfi = mean(median_mfi, na.rm = TRUE),
      sd_mfi = sd(median_mfi, na.rm = TRUE),
      cv_pct = 100 * sd_mfi / mean_mfi,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(cv_pct))
}

#' Bar chart of pool reproducibility (CV%) by antigen: one bar per antigen,
#' height = CV%, with a dashed reference line (default 20%) and bars above
#' the threshold colored red.
plot_pool_cv <- function(merged, pool_type = "pool", cv_threshold = 20,
                         exclude_controls = TRUE, controls = CONTROL_LABELS,
                         interactive = TRUE) {
  cv <- compute_pool_cv(merged, pool_type, exclude_controls, controls) |>
    dplyr::mutate(
      # .na_rm = FALSE: with few pool replicates (e.g. only 1 pool well on a
      # plate) sd()/cv_pct is NA for some or all antigens. The default
      # fct_reorder() drops NAs internally and errors out ("idx must contain
      # one integer for each level") when *every* level ends up NA-only;
      # keeping NAs just falls back to the antigens' original order instead
      # of crashing the plot.
      antigen = forcats::fct_reorder(antigen, cv_pct, .na_rm = FALSE),
      above = cv_pct > cv_threshold,
      tooltip = sprintf("%s\nCV: %.1f%%\nn = %d pools\nmean MFI: %.0f",
                        antigen, cv_pct, n, mean_mfi)
    )

  p <- suppressWarnings(
    ggplot(cv, aes(x = antigen, y = cv_pct)) +
      geom_col(aes(fill = above, text = tooltip)) +
      geom_hline(yintercept = cv_threshold, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "#2C7FB8"), guide = "none") +
      labs(x = NULL, y = "CV % (Median MFI, pools)",
           title = "Pool reproducibility (CV%) by antigen",
           subtitle = paste0("Dashed line = ", cv_threshold, "% threshold; red = above threshold")) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  if (interactive) to_plotly(p) else p
}


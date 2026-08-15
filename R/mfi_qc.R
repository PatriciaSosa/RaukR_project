# Median MFI QC: control-anchored thresholds, abundance plots, and pool CV%.
# Pure functions (no Shiny); interactive plots are ggiraph girafe objects.

library(tidyverse)
library(ggiraph)

# Default control-bead names, matching the assay design in the project plan.
# Override these if a real antigen table uses different labels.
CONTROL_LABELS <- c(empty = "Bare", pos_igg = "ahIgG",
                    neg_igm = "ahIgM", semi_pos = "EBNA1")

is_control_antigen <- function(antigen, controls = CONTROL_LABELS) {
  antigen %in% controls
}

# ---------------------------------------------------------------------------
# Control-anchored MFI thresholds
# ---------------------------------------------------------------------------

#' Minimum/maximum Median MFI thresholds anchored to control beads.
#'
#' min = mean MFI of the empty/Bare bead (background, expected in every well).
#' max = mean MFI of the ahIgG bead in samples/pools only (blanks have no real
#' antibody, so including them would understate the assay's achievable ceiling).
mfi_control_thresholds <- function(merged,
                                   empty_antigen = CONTROL_LABELS[["empty"]],
                                   pos_antigen = CONTROL_LABELS[["pos_igg"]],
                                   neg_sample_types = c("sample", "pool", "blank"),
                                   pos_sample_types = c("sample", "pool")) {
  min_thr <- merged |>
    dplyr::filter(antigen == empty_antigen, sample_type %in% neg_sample_types) |>
    dplyr::summarise(m = mean(median_mfi, na.rm = TRUE)) |>
    dplyr::pull(m)
  max_thr <- merged |>
    dplyr::filter(antigen == pos_antigen, sample_type %in% pos_sample_types) |>
    dplyr::summarise(m = mean(median_mfi, na.rm = TRUE)) |>
    dplyr::pull(m)
  list(min = min_thr, max = max_thr,
       empty_antigen = empty_antigen, pos_antigen = pos_antigen)
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

  p <- ggplot(df, aes(x = antigen, y = median_mfi, color = sample_type)) +
    (if (interactive) {
      geom_jitter_interactive(aes(tooltip = tooltip, data_id = well_384),
                              width = 0.15, alpha = 0.6, size = 1.4)
    } else geom_jitter(width = 0.15, alpha = 0.6, size = 1.4)) +
    scale_color_manual(values = SAMPLE_TYPE_COLORS, name = "Sample type") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Antigens", y = "Raw MFI Values",
         title = "Raw MFI value vs control beads") +
    theme_minimal(base_size = 11)

  if (interactive) girafe(ggobj = p, width_svg = 8, height_svg = 5) else p
}

# ---------------------------------------------------------------------------
# Boxplot: Median MFI by sample (identify low/high samples by clinical code)
# ---------------------------------------------------------------------------

#' Boxplot of Median MFI per sample (one box per sample_id across its
#' antigens), colored by sample type. Hovering a point reveals the clinical
#' code. Optionally restrict to a Quadrant/sample_type subset for readability.
plot_mfi_by_sample <- function(merged, quadrant = NULL, sample_types = NULL,
                               exclude_controls = TRUE, controls = CONTROL_LABELS,
                               thresholds = NULL, interactive = TRUE) {
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

  p <- ggplot(df, aes(x = sample_id, y = median_mfi)) +
    geom_boxplot(aes(color = sample_type), outlier.shape = NA, fill = NA, linewidth = 0.3) +
    (if (interactive) {
      geom_jitter_interactive(aes(color = sample_type, tooltip = tooltip, data_id = sample_id),
                              width = 0.2, alpha = 0.5, size = 0.9)
    } else geom_jitter(aes(color = sample_type), width = 0.2, alpha = 0.5, size = 0.9)) +
    scale_color_manual(values = SAMPLE_TYPE_COLORS, name = "Sample type") +
    scale_y_log10()

  if (!is.null(thresholds)) {
    p <- p +
      geom_hline(yintercept = thresholds$min, linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = thresholds$max, linetype = "dashed", color = "grey40")
  }

  p <- p +
    labs(x = "Cohort", y = "Signal (log scale)",
         title = "Raw MFI value vs all samples",
         subtitle = "Signals across samples") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (interactive) girafe(ggobj = p, width_svg = 11, height_svg = 5) else p
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

  p <- ggplot(df, aes(x = antigen, y = median_mfi)) +
    geom_boxplot(aes(fill = Antigen_Group), outlier.shape = NA, alpha = 0.3, linewidth = 0.3) +
    (if (interactive) {
      geom_jitter_interactive(aes(color = sample_type, tooltip = tooltip, data_id = well_384),
                              width = 0.15, alpha = 0.5, size = 0.8)
    } else geom_jitter(aes(color = sample_type), width = 0.15, alpha = 0.5, size = 0.8)) +
    scale_color_manual(values = SAMPLE_TYPE_COLORS, name = "Sample type") +
    scale_y_continuous(labels = scales::comma) +
    guides(fill = guide_legend(title = "Antigen group")) +
    labs(x = "Antigens", y = "Raw MFI Values",
         title = "Raw MFI value vs antigens") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (interactive) girafe(ggobj = p, width_svg = 11, height_svg = 5.5) else p
}

# ---------------------------------------------------------------------------
# Boxplot: Median MFI by sample_type, faceted by antigen (user-selectable)
# ---------------------------------------------------------------------------

plot_mfi_by_antigen_facet <- function(merged, antigens = NULL, interactive = TRUE) {
  df <- merged
  if (!is.null(antigens)) df <- dplyr::filter(df, antigen %in% antigens)
  df <- df |>
    dplyr::mutate(tooltip = sprintf("Sample: %s\nMFI: %.0f", sample_id, median_mfi))

  p <- ggplot(df, aes(x = sample_type, y = median_mfi, color = sample_type)) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    (if (interactive) {
      geom_jitter_interactive(aes(tooltip = tooltip, data_id = well_384),
                              width = 0.2, alpha = 0.4, size = 0.8)
    } else geom_jitter(width = 0.2, alpha = 0.4, size = 0.8)) +
    scale_color_manual(values = SAMPLE_TYPE_COLORS, guide = "none") +
    scale_y_log10() +
    facet_wrap(~antigen, scales = "free_y") +
    labs(x = NULL, y = "Median MFI (log scale)",
         title = "Median MFI by sample type, per antigen") +
    theme_minimal(base_size = 11)

  if (interactive) girafe(ggobj = p, width_svg = 10, height_svg = 6) else p
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
      antigen = forcats::fct_reorder(antigen, cv_pct),
      above = cv_pct > cv_threshold,
      tooltip = sprintf("%s\nCV: %.1f%%\nn = %d pools\nmean MFI: %.0f",
                        antigen, cv_pct, n, mean_mfi)
    )

  p <- ggplot(cv, aes(x = antigen, y = cv_pct)) +
    (if (interactive) {
      geom_col_interactive(aes(fill = above, tooltip = tooltip, data_id = antigen))
    } else geom_col(aes(fill = above))) +
    geom_hline(yintercept = cv_threshold, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "#2C7FB8"), guide = "none") +
    labs(x = NULL, y = "CV % (Median MFI, pools)",
         title = "Pool reproducibility (CV%) by antigen",
         subtitle = paste0("Dashed line = ", cv_threshold, "% threshold; red = above threshold")) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (interactive) girafe(ggobj = p, width_svg = 9, height_svg = 5) else p
}


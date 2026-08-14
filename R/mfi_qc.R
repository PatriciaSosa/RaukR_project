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
# Abundance vs. control-bead MFI scatter
# ---------------------------------------------------------------------------

#' Long data pairing each non-control analyte's abundance in a well with the
#' control-bead MFI values measured in that same well.
mfi_vs_controls_data <- function(merged,
                                 controls = CONTROL_LABELS[c("empty", "pos_igg", "neg_igm")]) {
  ctrl_wide <- merged |>
    dplyr::filter(antigen %in% controls) |>
    dplyr::select(well_384, antigen, median_mfi) |>
    tidyr::pivot_wider(names_from = antigen, values_from = median_mfi,
                       names_prefix = "ctrl_")

  merged |>
    dplyr::filter(!antigen %in% CONTROL_LABELS) |>
    dplyr::left_join(ctrl_wide, by = "well_384") |>
    tidyr::pivot_longer(dplyr::starts_with("ctrl_"),
                        names_to = "control_bead", values_to = "control_mfi") |>
    dplyr::mutate(control_bead = stringr::str_remove(control_bead, "^ctrl_"))
}

#' Scatter: antigen abundance (y) vs. control-bead MFI (x), one facet per
#' control bead, colored by sample type.
plot_mfi_vs_controls <- function(merged,
                                 controls = CONTROL_LABELS[c("empty", "pos_igg", "neg_igm")],
                                 interactive = TRUE) {
  df <- mfi_vs_controls_data(merged, controls) |>
    dplyr::mutate(tooltip = sprintf(
      "Sample: %s (%s)\nAntigen: %s\n%s MFI: %.0f\nAbundance: %.0f",
      sample_id, sample_type, antigen, control_bead, control_mfi, median_mfi
    ))

  p <- ggplot(df, aes(x = control_mfi, y = median_mfi, color = sample_type)) +
    (if (interactive) {
      geom_point_interactive(aes(tooltip = tooltip, data_id = paste(well_384, antigen)),
                             alpha = 0.5, size = 1.2)
    } else geom_point(alpha = 0.5, size = 1.2)) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = SAMPLE_TYPE_COLORS, name = "Sample type") +
    facet_wrap(~control_bead, scales = "free_x") +
    labs(x = "Control bead Median MFI (log scale)",
         y = "Antigen Median MFI (log scale)",
         title = "Antigen abundance vs. control bead MFI") +
    theme_minimal(base_size = 11)

  if (interactive) girafe(ggobj = p, width_svg = 10, height_svg = 4.3) else p
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
    labs(x = NULL, y = "Median MFI (log scale)",
         title = "Median MFI by sample") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))

  if (interactive) girafe(ggobj = p, width_svg = 11, height_svg = 5) else p
}

# ---------------------------------------------------------------------------
# Boxplot: Median MFI by antigen (samples overlaid, colored by sample type)
# ---------------------------------------------------------------------------

plot_mfi_by_antigen <- function(merged, exclude_controls = FALSE,
                                controls = CONTROL_LABELS, interactive = TRUE) {
  df <- merged
  if (exclude_controls) df <- dplyr::filter(df, !antigen %in% controls)
  df <- df |>
    dplyr::mutate(tooltip = sprintf("Sample: %s (%s)\nAntigen: %s\nMFI: %.0f",
                                    sample_id, sample_type, antigen, median_mfi))

  p <- ggplot(df, aes(x = antigen, y = median_mfi)) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "grey40", linewidth = 0.3) +
    (if (interactive) {
      geom_jitter_interactive(aes(color = sample_type, tooltip = tooltip, data_id = well_384),
                              width = 0.15, alpha = 0.4, size = 0.8)
    } else geom_jitter(aes(color = sample_type), width = 0.15, alpha = 0.4, size = 0.8)) +
    scale_color_manual(values = SAMPLE_TYPE_COLORS, name = "Sample type") +
    scale_y_log10() +
    labs(x = NULL, y = "Median MFI (log scale)",
         title = "Median MFI by antigen") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (interactive) girafe(ggobj = p, width_svg = 10, height_svg = 5) else p
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

#' CV% vs antigen for pool samples, with a reference line (default 20%).
#' Control beads are excluded by default: their CV is naturally inflated by
#' near-background noise and doesn't reflect antigen reproducibility.
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
    geom_hline(yintercept = cv_threshold, linetype = "dashed", color = "grey40") +
    (if (interactive) {
      geom_point_interactive(aes(color = above, tooltip = tooltip, data_id = antigen), size = 3)
    } else geom_point(aes(color = above), size = 3)) +
    scale_color_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "#2C7FB8"), guide = "none") +
    labs(x = NULL, y = "CV % (Median MFI, pools)",
         title = "Pool reproducibility (CV%) by antigen",
         subtitle = paste0("Dashed line = ", cv_threshold, "% threshold; red = above threshold")) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (interactive) girafe(ggobj = p, width_svg = 9, height_svg = 5) else p
}

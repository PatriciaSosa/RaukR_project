# 384-well plate visualization for bead-count QC.
# Pure functions (no Shiny); the interactive plot is a ggiraph girafe object.

library(tidyverse)
library(ggiraph)

# Quality ranges (from the plan):
#   Bad      bead count <= 35
#   Moderate 35 < bead count < 50
#   Good     bead count >= 50
QC_BAD_MAX  <- 35
QC_GOOD_MIN <- 50
QUALITY_LEVELS <- c("Bad", "Moderate", "Good")
QUALITY_COLORS <- c(Bad = "#D7301F", Moderate = "#FDBB84", Good = "#31A354")

# Aggregate views collapse the 16 beads in a well to one value.
AGG_VIEWS <- c("min", "median", "mean")

#' Classify bead counts into Bad / Moderate / Good.
classify_bead_count <- function(count, bad_max = QC_BAD_MAX, good_min = QC_GOOD_MIN) {
  factor(
    dplyr::case_when(
      is.na(count)        ~ NA_character_,
      count <= bad_max    ~ "Bad",
      count >= good_min   ~ "Good",
      TRUE                ~ "Moderate"
    ),
    levels = QUALITY_LEVELS
  )
}

#' Full 16x24 grid of 384-plate positions with numeric row/col coordinates.
plate_full_grid <- function() {
  tidyr::expand_grid(row = 1:16, col = 1:24) |>
    dplyr::mutate(
      row_letter = LETTERS[row],
      well_384 = paste0(row_letter, col)
    )
}

#' Summarise bead count to one value per well.
#'
#' Each well has one bead count per analyte. `view` selects how to collapse the
#' 16 beads to a single value per well:
#'   - "min"    the well's worst bead (default; conservative QC)
#'   - "median" the well's typical bead
#'   - "mean"   the well's average bead
#'   - an analyte column name (e.g. "Analyte.90") to show just that bead.
#' `limiting_analyte` records which bead is lowest in the well.
well_bead_summary <- function(merged, view = "min",
                              bad_max = QC_BAD_MAX, good_min = QC_GOOD_MIN) {
  if (view %in% AGG_VIEWS) {
    aggfun <- switch(view, min = min, median = stats::median, mean = mean)
    summ <- merged |>
      dplyr::group_by(well_384, sample_id, sample_type, Quadrant) |>
      dplyr::summarise(
        bead_count = aggfun(bead_count),
        limiting_analyte = .data$analyte[which.min(bead_count)],
        .groups = "drop"
      )
  } else {
    if (!view %in% merged$analyte) {
      stop("Unknown view '", view, "'. Use one of ",
           paste(AGG_VIEWS, collapse = "/"), " or a valid analyte name.")
    }
    summ <- merged |>
      dplyr::filter(.data$analyte == view) |>
      dplyr::group_by(well_384, sample_id, sample_type, Quadrant) |>
      dplyr::summarise(
        bead_count = dplyr::first(bead_count),
        limiting_analyte = dplyr::first(.data$analyte),
        .groups = "drop"
      )
  }

  summ |>
    dplyr::mutate(quality = classify_bead_count(bead_count, bad_max, good_min))
}

#' Named choices for a plate-view dropdown: the three aggregates plus one entry
#' per analyte labeled with its antigen name.
plate_view_choices <- function(merged) {
  ants <- merged |>
    dplyr::distinct(analyte, antigen) |>
    dplyr::arrange(antigen)
  analyte_choices <- stats::setNames(
    ants$analyte,
    ifelse(is.na(ants$antigen), ants$analyte,
           paste0(ants$antigen, " (", ants$analyte, ")"))
  )
  c(
    "Worst bead \u2014 min" = "min",
    "Median across beads"   = "median",
    "Mean across beads"     = "mean",
    analyte_choices
  )
}

#' Interactive 384-well plate colored by bead-count quality.
#'
#' Returns a ggiraph girafe object (hover tooltip + multiple click selection)
#' when `interactive = TRUE`, otherwise a plain ggplot (useful for exports).
plot_plate_beadcount <- function(merged, view = "min",
                                 bad_max = QC_BAD_MAX, good_min = QC_GOOD_MIN,
                                 interactive = TRUE) {
  summ <- well_bead_summary(merged, view, bad_max, good_min)
  is_agg <- view %in% AGG_VIEWS

  grid <- plate_full_grid() |>
    dplyr::left_join(summ, by = "well_384") |>
    dplyr::mutate(
      tooltip = dplyr::if_else(
        is.na(bead_count),
        paste0("Well ", well_384, "\n(no data)"),
        sprintf(
          "Well %s\nSample: %s (%s)\nBead count: %s\nQuality: %s%s",
          well_384, sample_id, sample_type, round(bead_count), quality,
          if (is_agg) paste0("\nLimiting bead: ", limiting_analyte) else ""
        )
      )
    )

  subtitle <- if (is_agg) {
    paste0("Colored by well ", view, " bead count across all beads")
  } else {
    ant_name <- merged$antigen[merged$analyte == view][1]
    label <- if (is.na(ant_name)) view else paste0(ant_name, " (", view, ")")
    paste0("Colored by bead count for ", label)
  }

  p <- ggplot(grid, aes(x = col, y = row)) +
    geom_tile_interactive(
      aes(fill = quality, tooltip = tooltip, data_id = well_384),
      color = "white", linewidth = 0.4
    ) +
    scale_fill_manual(
      values = QUALITY_COLORS, na.value = "grey90",
      drop = FALSE, name = "Bead count QC"
    ) +
    scale_x_continuous(breaks = 1:24, position = "top", expand = c(0, 0)) +
    scale_y_reverse(breaks = 1:16, labels = LETTERS[1:16], expand = c(0, 0)) +
    coord_equal() +
    labs(x = NULL, y = NULL,
         title = "384-well plate — bead count QC",
         subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  if (!interactive) return(p)

  girafe(
    ggobj = p, width_svg = 9, height_svg = 6.2,
    options = list(
      opts_hover(css = "stroke:black;stroke-width:1.5px;cursor:pointer;"),
      opts_selection(type = "multiple", css = "stroke:black;stroke-width:2px;"),
      opts_tooltip(css = "background:#333;color:#fff;padding:6px;border-radius:4px;font-size:12px;")
    )
  )
}

# Ingestion & parsing helpers for the Luminex FlexMap 3D serology QC app.
# Pure functions (no Shiny) so they can be unit-tested on their own.

library(tidyverse)
library(dplyr)
library(stats)
library(readr)
library(plotly)

# ---------------------------------------------------------------------------
# Raw xPONENT CSV parsing
# ---------------------------------------------------------------------------

# Split a single CSV line into fields, respecting double quotes.
.split_csv_line <- function(line) {
  scan(text = line, what = "character", sep = ",",
       quote = "\"", quiet = TRUE, strip.white = TRUE)
}

#' Parse a FlexMap 3D / xPONENT raw CSV export.
#'
#' Returns a list with:
#'   - metadata:  tibble(field, value) from the header block
#'   - tables:    named list of tibbles, one per DataType section
#'   - datatypes: character vector of DataType names found
#'
#' The export is a metadata header followed by `"DataType:","<name>"` sections,
#' each a wide table with `Location`, `Sample`, one column per analyte, and a
#' trailing `Total Events` column.
parse_flexmap_csv <- function(path) {
  lines <- readr::read_lines(path)
  blank <- trimws(lines) %in% c("", "\"\"")

  dt_lines <- which(str_detect(lines, '^\\s*"?DataType:"?'))
  if (length(dt_lines) == 0) {
    stop("No 'DataType:' sections found - is this a FlexMap 3D export?")
  }

  # ---- metadata: everything before the first DataType section ----
  meta_idx <- seq_len(dt_lines[1] - 1)
  meta_idx <- meta_idx[!blank[meta_idx]]
  metadata <- map_dfr(meta_idx, function(i) {
    f <- .split_csv_line(lines[i])
    if (length(f) == 0) return(NULL)
    tibble(field = f[1],
           value = if (length(f) > 1) paste(f[-1], collapse = ", ") else NA_character_)
  })

  # ---- each DataType section ----
  parse_section <- function(dt_start) {
    name <- .split_csv_line(lines[dt_start])[2]
    header_idx <- dt_start + 1
    # data rows run until the next blank line, next DataType, or CRC marker
    j <- header_idx + 1
    n <- length(lines)
    while (j <= n && !blank[j] &&
           !str_detect(lines[j], '^\\s*"?DataType:"?') &&
           !str_detect(lines[j], "CRC")) {
      j <- j + 1
    }
    block <- paste(lines[header_idx:(j - 1)], collapse = "\n")
    df <- readr::read_csv(I(block), show_col_types = FALSE,
                          name_repair = "minimal")
    list(name = name, data = df)
  }

  sections <- map(dt_lines, parse_section)
  tables <- set_names(map(sections, "data"), map_chr(sections, "name"))

  list(metadata = metadata,
       tables = tables,
       datatypes = names(tables))
}

#' Reshape the important Median + Count tables into one tidy long table.
#'
#' Location strings look like `"12(1,B3)"`: running index (plate, well_384).
#' Analyte column names look like `"Analyte.90"`, where the number is the
#' bead_id used to match the antigen table.
#' Output columns: run_idx, location, plate, well_384, raw_sample, analyte,
#' bead_id, median_mfi, bead_count, total_events.
flexmap_to_long <- function(parsed,
                            median_name = "Median",
                            count_name = "Count") {
  if (!median_name %in% names(parsed$tables)) {
    stop("Median table '", median_name, "' not found. Available: ",
         paste(parsed$datatypes, collapse = ", "))
  }
  if (!count_name %in% names(parsed$tables)) {
    stop("Count table '", count_name, "' not found. Available: ",
         paste(parsed$datatypes, collapse = ", "))
  }

  meta_cols <- c("Location", "Sample", "Total Events")

  tidy_one <- function(df, value_name) {
    analytes <- setdiff(names(df), meta_cols)
    df |>
      select(Location, Sample, `Total Events`, all_of(analytes)) |>
      pivot_longer(all_of(analytes), names_to = "analyte",
                   values_to = value_name)
  }

  med <- tidy_one(parsed$tables[[median_name]], "median_mfi")
  cnt <- tidy_one(parsed$tables[[count_name]], "bead_count") |>
    select(Location, analyte, bead_count)

  med |>
    left_join(cnt, by = c("Location", "analyte")) |>
    # Location "12(1,B3)" -> run_idx 12, plate 1, well_384 B3
    extract(Location, c("run_idx", "plate", "well_384"),
            regex = "^(\\d+)\\((\\d+),([A-P]\\d+)\\)$",
            remove = FALSE, convert = TRUE) |>
    # "Analyte.90" -> bead_id 90 (matches bead_id in the antigen table)
    mutate(bead_id = as.numeric(str_extract(analyte, "\\d+"))) |>
    rename(location = Location, raw_sample = Sample,
           total_events = `Total Events`) |>
    relocate(run_idx, location, plate, well_384, raw_sample, analyte,
             bead_id, median_mfi, bead_count, total_events)
}

# ---------------------------------------------------------------------------
# Templates the user downloads, completes, and re-uploads
# ---------------------------------------------------------------------------

ANTIGEN_COLS <- c("antigen", "bead_id", "virus_genus", "Pathogen", "Antigen_Group")
TRACE_COLS   <- c("well_96", "source_plate", "Quadrant", "well_384",
                  "sample_id", "sample_type")
SAMPLE_TYPES <- c("pool", "blank", "sample")

#' Empty antigen template pre-filled with the bead_ids seen in the raw data.
#' Bead ids are taken from the "Analyte.<bead_id>" column names; the user fills
#' in the antigen name and annotation for each.
antigen_template <- function(parsed = NULL) {
  bead_ids <- if (is.null(parsed)) numeric(0) else {
    df <- parsed$tables[[parsed$datatypes[1]]]
    cols <- setdiff(names(df), c("Location", "Sample", "Total Events"))
    # Real xPONENT/FlexMap exports name analyte columns "Analyte 255" (a
    # space), not "Analyte.255" (a dot) - read_csv() with name_repair =
    # "minimal" keeps the header exactly as exported. str_remove() with a
    # literal dot never matched real data, so bead_id came out as NA for
    # every row. Extract the digits directly instead, the same permissive
    # way flexmap_to_long() already does, so this works regardless of the
    # separator style a given export uses.
    as.numeric(str_extract(cols, "\\d+"))
  }
  tibble(bead_id = bead_ids) |>
    mutate(
      antigen = NA_character_,
      virus_genus = NA_character_,
      Pathogen = NA_character_,
      Antigen_Group = NA_character_
    ) |>
    select(all_of(ANTIGEN_COLS))
}

#' Empty traceability template pre-filled with the wells seen in the raw data.
traceability_template <- function(long = NULL) {
  wells <- if (is.null(long)) character(0) else sort(unique(long$well_384))
  tibble(
    well_96 = NA_character_,
    source_plate = NA_character_,
    Quadrant = NA_character_,
    well_384 = wells,
    sample_id = NA_character_,
    sample_type = NA_character_
  )
}

# ---------------------------------------------------------------------------
# Validation of uploaded tables
# ---------------------------------------------------------------------------

.new_report <- function() list(ok = TRUE, errors = character(0), warnings = character(0))
.err <- function(rep, msg) { rep$ok <- FALSE; rep$errors <- c(rep$errors, msg); rep }
.warn <- function(rep, msg) { rep$warnings <- c(rep$warnings, msg); rep }

#' Validate an uploaded antigen table. Returns list(ok, errors, warnings, data).
validate_antigen_table <- function(df, long = NULL) {
  rep <- .new_report()
  missing_cols <- setdiff(ANTIGEN_COLS, names(df))
  if (length(missing_cols)) {
    rep <- .err(rep, paste("Missing columns:", paste(missing_cols, collapse = ", ")))
    return(c(rep, list(data = df)))
  }

  bead <- suppressWarnings(as.numeric(df$bead_id))
  if (any(is.na(bead) & !is.na(df$bead_id)))
    rep <- .err(rep, "bead_id must be numeric.")
  if (any(bead < 1 | bead > 384, na.rm = TRUE))
    rep <- .err(rep, "bead_id must be between 1 and 384.")
  if (any(is.na(df$antigen)) )
    rep <- .err(rep, "antigen has missing values.")
  if (anyDuplicated(df$antigen))
    rep <- .err(rep, "antigen values must be unique.")
  if (anyDuplicated(na.omit(bead)))
    rep <- .warn(rep, "Duplicate bead_id values present.")

  if (!is.null(long)) {
    raw_beads <- unique(long$bead_id)
    unmapped <- setdiff(raw_beads, bead)
    if (length(unmapped))
      rep <- .warn(rep, paste("Bead IDs in raw data not in antigen table:",
                              paste(unmapped, collapse = ", ")))
  }
  c(rep, list(data = mutate(df, bead_id = bead)))
}

#' Validate an uploaded traceability table. Returns list(ok, errors, warnings, data).
validate_traceability_table <- function(df, long = NULL) {
  rep <- .new_report()
  missing_cols <- setdiff(TRACE_COLS, names(df))
  if (length(missing_cols)) {
    rep <- .err(rep, paste("Missing columns:", paste(missing_cols, collapse = ", ")))
    return(c(rep, list(data = df)))
  }

  # Real hand-edited CSVs commonly carry stray whitespace or inconsistent
  # case ("Pool", " blank "). merge_qc_data() later does
  # factor(sample_type, levels = SAMPLE_TYPES): any value that isn't an
  # *exact* match to "pool"/"blank"/"sample" is silently turned into NA,
  # which looks like "the app doesn't recognize pools and blanks" with no
  # error message anywhere. Normalize here, before validation, so this
  # class of formatting slip is either fixed automatically or caught as an
  # explicit error below instead of failing silently downstream.
  df <- df |>
    dplyr::mutate(
      sample_type = tolower(trimws(sample_type)),
      well_384 = toupper(trimws(well_384))
    )

  if (any(is.na(df$sample_type)))
    rep <- .warn(rep, "Some sample_type values are missing; those wells will not be classified as sample/pool/blank.")

  bad_type <- setdiff(unique(na.omit(df$sample_type)), SAMPLE_TYPES)
  if (length(bad_type))
    rep <- .err(rep, paste("Invalid sample_type values:", paste(bad_type, collapse = ", "),
                           "- allowed:", paste(SAMPLE_TYPES, collapse = ", ")))
  if (any(!str_detect(df$well_384, "^[A-P]([1-9]|1[0-9]|2[0-4])$"), na.rm = TRUE))
    rep <- .err(rep, "well_384 has values outside a 16x24 (A-P, 1-24) plate.")
  if (anyDuplicated(df$well_384))
    rep <- .err(rep, "well_384 values must be unique.")
  if (any(is.na(df$sample_id)))
    rep <- .warn(rep, "Some sample_id values are missing.")
  bad_quad <- setdiff(unique(na.omit(df$Quadrant)), paste0("Q", 1:4))
  if (length(bad_quad))
    rep <- .warn(rep, paste("Unexpected Quadrant values:", paste(bad_quad, collapse = ", ")))

  if (!is.null(long)) {
    raw_wells <- unique(long$well_384)
    unmapped <- setdiff(raw_wells, df$well_384)
    if (length(unmapped))
      rep <- .warn(rep, paste(length(unmapped),
                              "well(s) in raw data not in traceability table."))
  }
  c(rep, list(data = df))
}

# ---------------------------------------------------------------------------
# Merge everything into the analysis-ready long table
# ---------------------------------------------------------------------------

#' Join the tidy raw long table with the antigen and traceability annotations.
#'
#'   raw + antigen table:  Analyte.<bead_id>  <->  antigen table `bead_id`
#'   + traceability:       well_384 (from Location)  <->  traceability `well_384`
#'
#' Authoritative sample_id and sample_type come from the traceability table
#' (keyed by well_384); the raw Sample label is retained as raw_sample for
#' cross-checking.
merge_qc_data <- function(long, antigen_tbl, trace_tbl) {
  long |>
    left_join(antigen_tbl, by = "bead_id") |>
    left_join(trace_tbl, by = "well_384") |>
    mutate(sample_type = factor(sample_type, levels = SAMPLE_TYPES)) |>
    relocate(sample_id, sample_type, source_plate, Quadrant, well_96,
             well_384, antigen, bead_id, median_mfi, bead_count)
}

# Colors for sample_type, reused across plots.
SAMPLE_TYPE_COLORS <- c(sample = "#1239c4", pool = "#c55905", blank = "#393838")

#' Convert a ggplot to an interactive plotly widget with sensible defaults:
#' hover tooltips sourced from the `text` aesthetic, and the mode bar enabled
#' for zoom / pan / box-select / reset-axes (plotly's built-in interactions).
to_plotly <- function(p, tooltip = "text") {
  plotly::ggplotly(p, tooltip = tooltip) |>
    plotly::layout(dragmode = "zoom", hoverlabel = list(bgcolor = "white")) |>
    plotly::config(displaylogo = FALSE)
}

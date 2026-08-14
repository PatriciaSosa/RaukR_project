# Manual tests for R/ingest.R against the synthetic data.
library(tidyverse)
proj <- "/Users/patriciasosa/Documents/RaukR_proj/project"
source(file.path(proj, "R/ingest.R"))

pass <- function(label, cond) {
  cat(if (isTRUE(cond)) "PASS  " else "FAIL  ", label, "\n")
  invisible(cond)
}

# ---- 1. Parse raw CSV ----
parsed <- parse_flexmap_csv(file.path(proj, "data/raw/flexmap_raw.csv"))
pass("DataTypes found: Median, Net MFI, Count, Mean",
     setequal(parsed$datatypes, c("Median", "Net MFI", "Count", "Mean")))
pass("Instrument metadata captured",
     any(parsed$metadata$field == "Instrument Type"))
pass("Median table is 384 x 19",
     all(dim(parsed$tables$Median) == c(384, 19)))
pass("Raw analyte columns named Analyte.<bead_id>",
     all(str_detect(setdiff(names(parsed$tables$Median),
                            c("Location", "Sample", "Total Events")),
                    "^Analyte\\.\\d+$")))

# ---- 2. Tidy to long ----
long <- flexmap_to_long(parsed)
pass("Long table has 384*16 = 6144 rows", nrow(long) == 384 * 16)
pass("Location parsed into well_384", all(!is.na(long$well_384)))
pass("bead_id extracted from Analyte.X column names", all(!is.na(long$bead_id)))
pass("median_mfi and bead_count numeric",
     is.numeric(long$median_mfi) && is.numeric(long$bead_count))
pass("bead_count join complete (no NA)", !any(is.na(long$bead_count)))

# ---- 3. Templates ----
at <- antigen_template(parsed)
pass("Antigen template pre-fills 16 bead_ids", nrow(at) == 16 && all(!is.na(at$bead_id)))
pass("Antigen template has required cols", all(ANTIGEN_COLS %in% names(at)))
tt <- traceability_template(long)
pass("Traceability template pre-fills 384 wells", nrow(tt) == 384)
pass("Traceability template has required cols", all(TRACE_COLS %in% names(tt)))

# ---- 4. Validation of the (good) filled templates ----
ant  <- read_csv(file.path(proj, "data/templates/antigen_table.csv"), show_col_types = FALSE)
trace <- read_csv(file.path(proj, "data/templates/traceability_table.csv"), show_col_types = FALSE)
va <- validate_antigen_table(ant, long)
vt <- validate_traceability_table(trace, long)
pass("Good antigen table validates ok", va$ok && length(va$errors) == 0)
pass("Good traceability table validates ok", vt$ok && length(vt$errors) == 0)

# ---- 5. Validation catches broken inputs ----
bad_ant <- ant; bad_ant$bead_id[1] <- 999
pass("Out-of-range bead_id flagged", !validate_antigen_table(bad_ant)$ok)
missing_ant <- ant |> select(-Pathogen)
pass("Missing antigen column flagged", !validate_antigen_table(missing_ant)$ok)
bad_trace <- trace; bad_trace$sample_type[1] <- "positive"
pass("Invalid sample_type flagged", !validate_traceability_table(bad_trace)$ok)
dup_trace <- trace; dup_trace$well_384[2] <- dup_trace$well_384[1]
pass("Duplicate well_384 flagged", !validate_traceability_table(dup_trace)$ok)

# ---- 6. Merge ----
merged <- merge_qc_data(long, va$data, vt$data)
pass("Merged rows preserved (6144)", nrow(merged) == 6144)
pass("sample_id populated from traceability (Location->well_384)", !any(is.na(merged$sample_id)))
pass("antigen name populated from antigen-table join (Analyte.X->bead_id)", !any(is.na(merged$antigen)))
pass("bead_id present throughout", !any(is.na(merged$bead_id)))
pass("sample_type is a factor with 3 levels",
     is.factor(merged$sample_type) && nlevels(merged$sample_type) == 3)
pass("raw_sample matches traceability sample_id",
     all(merged$raw_sample == merged$sample_id))

cat("\n--- merged glimpse ---\n")
glimpse(merged)
cat("\n--- sample_type x quadrant ---\n")
print(distinct(merged, sample_id, sample_type, Quadrant) |> count(sample_type, Quadrant))

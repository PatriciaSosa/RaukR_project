# Manual tests for R/refine.R against the synthetic data.
library(tidyverse)
proj <- "/Users/patriciasosa/Documents/RaukR_proj/project"
source(file.path(proj, "R/ingest.R"))
source(file.path(proj, "R/plate_view.R"))
source(file.path(proj, "R/mfi_qc.R"))
source(file.path(proj, "R/refine.R"))

pass <- function(label, cond) {
  cat(if (isTRUE(cond)) "PASS  " else "FAIL  ", label, "\n")
  invisible(cond)
}

parsed <- parse_flexmap_csv(file.path(proj, "data/raw/flexmap_raw.csv"))
long <- flexmap_to_long(parsed)
ant  <- read_csv(file.path(proj, "data/templates/antigen_table.csv"), show_col_types = FALSE)
trace <- read_csv(file.path(proj, "data/templates/traceability_table.csv"), show_col_types = FALSE)
merged <- merge_qc_data(long, validate_antigen_table(ant)$data, validate_traceability_table(trace)$data)

# ---- 1. Flags table ----
flags <- build_sample_flags(merged)
pass("One flag row per well (384)", nrow(flags) == 384)
pass("Bead-count categories match Phase 3 (10 Bad, 15 Moderate, 359 Good)",
     all(count(flags, bead_quality)$n == c(10, 15, 359)))
cat("\nflag_reason distribution:\n"); print(count(flags, flag_reason))

# ---- 2. flagged_for_review ----
review <- flagged_for_review(flags)
pass("Review subset = Bad + Moderate = 25",
     nrow(review) == sum(flags$bead_quality != "Good"))

# ---- 3. Manual removal list round-trip ----
manual_ids <- flags |> filter(bead_quality == "Good") |> slice_sample(n = 3) |> pull(sample_id)
tmp <- tempfile(fileext = ".csv")
write_csv(tibble(sample_id = manual_ids), tmp)
read_back <- read_manual_removal_list(tmp)
pass("Manual removal list round-trips", setequal(read_back, manual_ids))

# ---- 4. Apply removals: bead count Bad only ----
res_bad <- apply_qc_removals(merged, flags, remove_bead_categories = "Bad")
pass("Bad-only removal drops exactly 10 wells", n_distinct(res_bad$removed$well_384) == 10)
pass("Clean table rows = (384-10)*16", nrow(res_bad$clean) == (384 - 10) * 16)
pass("No Bad wells remain in clean data",
     !any(res_bad$clean$well_384 %in% filter(flags, bead_quality == "Bad")$well_384))

# ---- 5. Apply removals: Bad + Moderate + manual ----
res_all <- apply_qc_removals(merged, flags,
                             remove_bead_categories = c("Bad", "Moderate"),
                             manual_sample_ids = manual_ids)
n_expected_dropped <- n_distinct(filter(flags, bead_quality %in% c("Bad", "Moderate"))$well_384) +
  length(manual_ids)  # manual ids are Good-quality, so no overlap
pass("Bad+Moderate+manual removal drops the right count",
     n_distinct(res_all$removed$well_384) == n_expected_dropped)
pass("Removed table carries a removal_reason for every row", !any(is.na(res_all$removed$removal_reason)))
pass("Log has 4 lines summarizing the run", length(res_all$log) == 4)
cat("\n--- removal log ---\n"); cat(paste(res_all$log, collapse = "\n"), "\n")

# ---- 6. No removals requested -> clean == merged ----
res_none <- apply_qc_removals(merged, flags)
pass("No removals leaves data untouched", nrow(res_none$clean) == nrow(merged))

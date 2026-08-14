# Generate synthetic Luminex FlexMap 3D (xPONENT-style) data for app development.
# Produces:
#   data/raw/flexmap_raw.csv        - raw equipment export (Median, Count, Net MFI, Mean)
#   data/templates/antigen_table.csv        - filled antigen annotation (example)
#   data/templates/traceability_table.csv   - filled traceability (example)
# The templates are also what the app would hand the user as blanks to fill in.

library(tidyverse)
library(glue)

set.seed(4173)

proj <- "/Users/patriciasosa/Documents/RaukR_proj/project"

# ---- Analyte / bead definitions -------------------------------------------
# Viral antigens + control beads. bead_id is a number in 1-384.
analytes <- tribble(
  ~antigen,      ~bead_id, ~virus_genus,     ~Pathogen,        ~Antigen_Group, ~role,
  "DENV1_NS1",   12,       "Flavivirus",     "Dengue virus 1", "NS1",          "antigen",
  "DENV2_NS1",   18,       "Flavivirus",     "Dengue virus 2", "NS1",          "antigen",
  "DENV3_NS1",   25,       "Flavivirus",     "Dengue virus 3", "NS1",          "antigen",
  "DENV4_NS1",   33,       "Flavivirus",     "Dengue virus 4", "NS1",          "antigen",
  "ZIKV_NS1",    42,       "Flavivirus",     "Zika virus",     "NS1",          "antigen",
  "ZIKV_EDIII",  45,       "Flavivirus",     "Zika virus",     "Envelope",     "antigen",
  "CHIKV_E1",    51,       "Alphavirus",     "Chikungunya",    "Envelope",     "antigen",
  "CHIKV_E2",    56,       "Alphavirus",     "Chikungunya",    "Envelope",     "antigen",
  "WNV_NS1",     62,       "Flavivirus",     "West Nile virus","NS1",          "antigen",
  "YFV_NS1",     70,       "Flavivirus",     "Yellow fever",   "NS1",          "antigen",
  "SARS2_S1",    77,       "Betacoronavirus","SARS-CoV-2",     "Spike",        "antigen",
  "SARS2_N",     85,       "Betacoronavirus","SARS-CoV-2",     "Nucleocapsid", "antigen",
  "EBNA1",       20,       "Lymphocryptovirus","Epstein-Barr", "Nuclear",      "semi_pos_ctrl",
  "Bare",        90,       "Control",        "None",           "Empty_bead",   "neg_ctrl_empty",
  "ahIgG",       95,       "Control",        "None",           "Anti-human_IgG","pos_ctrl",
  "ahIgM",       100,      "Control",        "None",           "Anti-human_IgM","neg_ctrl"
)

# ---- Plate layout: one 384 plate = 4 x 96 source plates -------------------
# Standard quadrant interleave for 96 -> 384:
#   source plate 1 -> odd 384 row, odd 384 col
#   source plate 2 -> odd 384 row, even 384 col
#   source plate 3 -> even 384 row, odd 384 col
#   source plate 4 -> even 384 row, even 384 col
rows96 <- LETTERS[1:8]
cols96 <- 1:12
rows384 <- LETTERS[1:16]
cols384 <- 1:24

make_quadrant <- function(src_plate) {
  row_off <- if (src_plate %in% c(1, 2)) 0 else 1
  col_off <- if (src_plate %in% c(1, 3)) 0 else 1
  expand_grid(r96 = seq_len(8), c96 = seq_len(12)) |>
    mutate(
      source_plate = glue("P{src_plate}"),
      # SELMA liquid handler combines each 96 plate into one 384 quadrant
      Quadrant = glue("Q{src_plate}"),
      well_96 = glue("{rows96[r96]}{c96}"),
      r384 = 2 * (r96 - 1) + 1 + row_off,
      c384 = 2 * (c96 - 1) + 1 + col_off,
      well_384 = glue("{rows384[r384]}{c384}")
    )
}

layout <- map(1:4, make_quadrant) |> list_rbind()

# ---- Assign sample types & clinical codes ---------------------------------
# Per source plate: A1 & A2 = blanks, B1 & B2 = pools, remaining = samples.
layout <- layout |>
  arrange(source_plate, r96, c96) |>
  group_by(source_plate) |>
  mutate(
    sample_type = case_when(
      well_96 %in% c("A1", "A2") ~ "blank",
      well_96 %in% c("B1", "B2") ~ "pool",
      TRUE ~ "sample"
    )
  ) |>
  ungroup()

# Clinical sample_id codes
layout <- layout |>
  mutate(
    sample_id = case_when(
      sample_type == "blank" ~ glue("BLANK_{source_plate}_{row_number()}"),
      sample_type == "pool"  ~ glue("POOL_{source_plate}_{row_number()}"),
      TRUE ~ NA_character_
    )
  )
n_samp <- sum(layout$sample_type == "sample")
layout$sample_id[layout$sample_type == "sample"] <-
  sprintf("CLIN-%04d", seq_len(n_samp))

# Running well order across the 384 plate (row-major) for the Location field
layout <- layout |>
  mutate(
    r384n = match(str_extract(well_384, "^[A-P]"), rows384),
    c384n = as.integer(str_extract(well_384, "\\d+$"))
  ) |>
  arrange(r384n, c384n) |>
  mutate(run_idx = row_number())

# ---- Simulate MFI per well x analyte --------------------------------------
# Baseline expectation of Median MFI depending on role and sample type.
sim_median <- function(role, sample_type) {
  base <- switch(role,
    pos_ctrl = if (sample_type == "blank") rlnorm(1, log(120), 0.3)
               else rlnorm(1, log(22000), 0.15),
    neg_ctrl_empty = rlnorm(1, log(35), 0.4),
    neg_ctrl = rlnorm(1, log(70), 0.4),
    semi_pos_ctrl = if (sample_type == "blank") rlnorm(1, log(60), 0.4)
                    else rlnorm(1, log(3500), 0.6),
    antigen = if (sample_type == "blank") rlnorm(1, log(45), 0.4)
              else if (sample_type == "pool") rlnorm(1, log(2500), 0.25)
              else {
                # samples: mixture of seronegative (low) and seropositive (high)
                if (runif(1) < 0.55) rlnorm(1, log(80), 0.5)
                else rlnorm(1, log(4000), 0.8)
              }
  )
  round(base, 1)
}

# Bead count per well x analyte: mostly good, some moderate/bad wells.
# Give a handful of wells a systematically low bead count (aspiration issues).
bad_wells <- sample(layout$well_384, 10)
mod_wells <- sample(setdiff(layout$well_384, bad_wells), 15)

sim_count <- function(well_384) {
  if (well_384 %in% bad_wells) round(runif(1, 5, 34))
  else if (well_384 %in% mod_wells) round(runif(1, 36, 49))
  else round(rnorm(1, 90, 20)) |> max(50)
}

grid <- layout |>
  select(run_idx, well_384, sample_id, sample_type, source_plate) |>
  cross_join(analytes |> select(antigen, role))

grid <- grid |>
  rowwise() |>
  mutate(
    Median = sim_median(role, sample_type),
    Count  = sim_count(well_384)
  ) |>
  ungroup() |>
  mutate(
    NetMFI = round(pmax(Median - rlnorm(n(), log(30), 0.3), 0), 1),
    Mean   = round(Median * rlnorm(n(), log(1.05), 0.05), 1)
  )

# ---- Build wide tables per DataType ---------------------------------------
# Real FlexMap exports label bead columns "Analyte.<bead_id>", not antigen names.
analyte_cols <- paste0("Analyte.", analytes$bead_id)

to_wide <- function(value_col) {
  grid |>
    mutate(
      Location = sprintf("%d(1,%s)", run_idx, well_384),
      Sample = sample_id
    ) |>
    left_join(analytes |> select(antigen, bead_id), by = "antigen") |>
    mutate(analyte = paste0("Analyte.", bead_id)) |>
    select(run_idx, Location, Sample, analyte, val = all_of(value_col)) |>
    pivot_wider(names_from = analyte, values_from = val) |>
    arrange(run_idx) |>
    select(Location, Sample, all_of(analyte_cols)) |>
    mutate(`Total Events` = round(runif(n(), 800, 1600)))
}

wide_median <- to_wide("Median")
wide_count  <- to_wide("Count")
wide_net    <- to_wide("NetMFI")
wide_mean   <- to_wide("Mean")

# ---- Write raw csv in xPONENT-like layout ---------------------------------
raw_path <- file.path(proj, "data/raw/flexmap_raw.csv")
con <- file(raw_path, open = "w")
wl <- function(...) writeLines(paste(..., sep = ""), con)

# Metadata header block
wl('"Program","xPONENT"')
wl('"Build","4.3.229.0"')
wl('"Date","08/14/2026","Time","14:32:11"')
wl('"SN","FM3D12345"')
wl('"Session","1"')
wl('"Operator","p.sosa"')
wl('"TemplateID",""')
wl('"TemplateName","SerologyQC_16plex"')
wl('"SampleVolume","50 uL"')
wl('"BatchAuthor","p.sosa"')
wl('"BatchDescription","Viral surveillance serology panel - synthetic"')
wl('"ProtocolName","IgG_16plex"')
wl('"ProtocolVersion","1"')
wl('"ProtocolPlate","Plate","1","1"')
wl('"Instrument Type","FLEXMAP 3D"')
wl('""')
wl('"Samples","', nrow(wide_median), '","Min Events","0"')
wl('""')

write_section <- function(dtype, df) {
  wl('"DataType:","', dtype, '"')
  # header row
  writeLines(paste0('"', paste(names(df), collapse = '","'), '"'), con)
  # data rows
  for (i in seq_len(nrow(df))) {
    vals <- df[i, ]
    loc <- paste0('"', vals[["Location"]], '"')
    smp <- paste0('"', vals[["Sample"]], '"')
    nums <- vals[setdiff(names(df), c("Location", "Sample"))] |>
      unlist() |> as.character()
    writeLines(paste(c(loc, smp, nums), collapse = ","), con)
  }
  wl('""')
}

write_section("Median", wide_median)
write_section("Net MFI", wide_net)
write_section("Count", wide_count)
write_section("Mean", wide_mean)
wl('"-- CRC --"')
close(con)
message("Wrote raw csv: ", raw_path)

# ---- Write filled template tables -----------------------------------------
antigen_table <- analytes |>
  select(antigen, bead_id, virus_genus, Pathogen, Antigen_Group)
write_csv(antigen_table, file.path(proj, "data/templates/antigen_table.csv"))

traceability_table <- layout |>
  arrange(source_plate, well_96) |>
  select(well_96, source_plate, Quadrant, well_384, sample_id, sample_type)
write_csv(traceability_table, file.path(proj, "data/templates/traceability_table.csv"))

message("Wrote templates: antigen_table.csv, traceability_table.csv")
message("Bad-count wells: ", paste(sort(bad_wells), collapse = ", "))
message("Moderate-count wells: ", paste(sort(mod_wells), collapse = ", "))

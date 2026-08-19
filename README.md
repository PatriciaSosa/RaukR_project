# Serology Quality Control Shiny App for Affinity Proteomics analysis

This is a Shiny app for quality control (QC) of multiplexed bead-based serology
assays run on a **Luminex FlexMap 3D** instrument. The assay couples viral
antigens to color-coded magnetic beads and incubates them with plasma/serum
samples to capture virus-specific IgG antibodies to study viral surveillance.

The app evaluates, filters, and validates the raw FlexMap 3D export.
It assesses samples, pools, and blanks using **bead count** and **median MFI
signal**, visualizes the 384-well plate, and lets the analyst remove
low-quality samples before exporting a clean, tidy long format dataset and a QC
report documenting every decision as well as analytical specifications.

## Assay design

- **Sample types:** `sample` (biological plasma), `pool` (commercial plasma for
  assay performance), `blank` (bead pool incubated with PBS).
- **Control beads:** `Bare` (empty, negative control), `ahIgG` (positive
  control), `ahIgM` (negative control), `EBNA1` (semi-positive control).
- **Plate format:** four 96-well source plates are combined by the **SELMA
  liquid handler** into one 384-well Greiner plate (16 rows A–P × 24 columns),
  one quadrant (Q1–Q4) per source plate.

## Inputs

### 1. Raw Luminex export (CSV)
The FlexMap 3D / xPONENT export in wide format: a metadata header (instrument,
batch, protocol) followed by `DataType:` sections. The app uses the **Median**
(signal) and **Count** (bead count) tables; other tables (Net MFI, Mean) are
ignored.

### 2. Tables the user provides / completes
After the raw CSV is loaded, the app offers two blank templates to download,
fill in, and re-upload:

| Table | Required columns | Purpose |
|-------|-----------------|---------|
| **Antigen table** | `antigen`, `bead_id` (1–384), `virus_genus`, `Pathogen`, `Antigen_Group` | Maps each analyte (bead) to its antigen and annotation |
| **Traceability table** | `well_96`, `source_plate`, `Quadrant` (Q1–Q4), `well_384`, `sample_id` (clinical code), `sample_type` (pool/blank/sample) | Links each 384 well to its source plate/quadrant, clinical code, and sample type |

Filled example versions live in [`data/templates/`](data/templates/).

## QC workflow

1. Upload raw CSV → map `Location`/`Sample` and analyte columns.
2. Match analytes to the antigen table (bead IDs).
3. Match wells to the traceability table (clinical codes, quadrants).
4. Define controls (blanks, pools) and set QC thresholds.
5. Assess and visualize the data, review flags manually.
6. Remove poor samples.
7. Export the cleaned long-format table and QC report.

### Bead-count evaluation
Interactive 384-well plate view, wells colored by quality range:

| Category | Bead count |
|----------|-----------|
| Bad | ≤ 35 |
| Moderate | 35 < count < 50 |
| Good | ≥ 50 |

The analyst can hover to inspect wells and choose to drop Bad and/or Moderate
wells, or keep everything.

### Median MFI evaluation
Interactive plots include:

- Scatter of abundance vs. control-bead MFI (empty, ahIgG, ahIgM), colored by
  sample type.
- Boxplot of Median MFI by sample, colored by sample type.
- Boxplot of Median MFI by antigen with observations, colored by sample type.
- Boxplot of Median MFI by sample type, faceted/filterable by antigen.
- CV% vs. antigen for **pools only** (assay reproducibility), with a 20%
  reference line.

Analysts may also upload a list of `sample_id`s to remove (e.g., lab errors, outliers).

## Outputs

- **Cleaned dataset** in long format (`sample_id`, `sample_type`, `antigen`,
  `well_384`, `Quadrant`, Median MFI, bead count, annotations).
- **QC flags** — samples with bad or moderate bead count in at least one bead ID.
- **QC report** — figures, thresholds, removals, and equipment metadata from the
  raw export, documenting every operation applied to the data.
- **Figures** — downloadable QC plots.

## Project structure

```
project/
├── data-raw/generate_synthetic_data.R   # reproducible synthetic-data generator
├── data/
│   ├── raw/flexmap_raw.csv               # example FlexMap 3D export (384 wells, 16 analytes)
│   └── templates/                        # filled antigen + traceability examples
├── R/                                    # app modules and helper functions
├── exports/                              # generated outputs
└── www/                                  # static assets
```

## Example figures

![Plate view](www/Plate_View.png)
![Pool reproducibility](www/Pool_CV.png)

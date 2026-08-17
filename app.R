# Quality Control Shiny App for Affinity Proteomics - Serology Data
# Run from the project root: shiny::runApp(".")

library(shiny)
library(bslib)
library(tidyverse)
library(ggiraph)
library(DT)

source("R/ingest.R")
source("R/plate_view.R")
source("R/mfi_qc.R")
source("R/refine.R")
source("R/export.R")

DEFAULT_RAW   <- "data/raw/flexmap_raw.csv"
DEFAULT_AG    <- "data/templates/antigen_table.csv"
DEFAULT_TRACE <- "data/templates/traceability_table.csv"

navbar_css <- tags$style(HTML("
  .navbar { background-color: #b7e4b7 !important; }
  .navbar .navbar-brand, .navbar .nav-link { color: #000000 !important; font-weight: 700 !important; }
  .navbar .nav-link.active { text-decoration: underline; }
"))

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- tagList(
  tags$head(navbar_css),
  page_navbar(
    title = "Serology QC \u2014 Affinity Proteomics",
    theme = bs_theme(version = 5, primary = "#2C7FB8"),
    sidebar = sidebar(
      title = "1. Load data",
      fileInput("raw_csv", "Raw FlexMap 3D CSV", accept = ".csv"),
      hr(),
      downloadButton("dl_antigen_template", "Antigen template", class = "btn-sm w-100 mb-1"),
      downloadButton("dl_trace_template", "Traceability template", class = "btn-sm w-100 mb-2"),
      fileInput("antigen_csv", "Completed antigen table", accept = ".csv"),
      fileInput("trace_csv", "Completed traceability table", accept = ".csv"),
      hr(),
      uiOutput("data_status"),
      helpText("Until you upload your own files, the app runs on the bundled synthetic example data.")
    ),
    nav_panel(
      "Bead count QC",
      card(
        full_screen = TRUE,
        card_header("384-well plate \u2014 bead count"),
        selectInput("plate_view_mode", "View", choices = c("min"), selected = "min", width = "320px"),
        girafeOutput("plate_plot", height = "650px")
      )
    ),
    nav_panel(
      "MFI QC",
      layout_column_wrap(
        width = 1/2,
        card(full_screen = TRUE, card_header("Raw MFI value vs control beads"),
             girafeOutput("mfi_controls", height = "480px")),
        card(full_screen = TRUE, card_header("Raw MFI value vs antigens"),
             girafeOutput("mfi_antigen", height = "480px")),
        card(full_screen = TRUE, card_header("Raw MFI value vs all samples"),
             girafeOutput("mfi_sample", height = "480px")),
        card(full_screen = TRUE, card_header("Pool reproducibility (CV%)"),
             girafeOutput("mfi_cv", height = "480px"))
      ),
      card(
        full_screen = TRUE,
        card_header("Median MFI by sample type, for one antigen"),
        selectInput("facet_antigen", "Antigen", choices = NULL, width = "320px"),
        girafeOutput("mfi_facet", height = "550px")
      )
    ),
    nav_panel(
      "Refinement",
      layout_column_wrap(
        width = 1/2,
        card(
          card_header("QC flags (review before removing)"),
          DTOutput("flags_table")
        ),
        card(
          card_header("Removal decisions"),
          checkboxGroupInput("remove_categories", "Remove wells with bead-count quality:",
                             choices = c("Bad", "Moderate"), selected = character(0)),
          p(class = "text-muted small",
            "MFI-based removal is always the analyst's decision, based on reviewing the MFI plots \u2014 there is no automatic threshold removal here."),
          downloadButton("dl_sample_removal_template", "Download sample list template",
                         class = "btn-sm w-100 mb-2"),
          helpText("Edit the downloaded file down to just the sample_id rows you want removed, then upload it below."),
          fileInput("manual_removal_csv", "Manual removal list (CSV with sample_id column)", accept = ".csv"),
          actionButton("apply_removal", "Apply removal", class = "btn-primary"),
          hr(),
          verbatimTextOutput("removal_log")
        )
      )
    ),
    nav_panel(
      "Export",
      card(
        card_header("Downloads"),
        downloadButton("dl_clean", "Download cleaned dataset (long format)", class = "w-100 mb-2"),
        downloadButton("dl_flags", "Download QC flags", class = "w-100 mb-2"),
        downloadButton("dl_report", "Download QC report (PDF)", class = "w-100")
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  loaded <- reactive({
    raw_path <- if (!is.null(input$raw_csv)) input$raw_csv$datapath else DEFAULT_RAW
    ag_path  <- if (!is.null(input$antigen_csv)) input$antigen_csv$datapath else DEFAULT_AG
    tr_path  <- if (!is.null(input$trace_csv)) input$trace_csv$datapath else DEFAULT_TRACE

    parsed <- parse_flexmap_csv(raw_path)
    long <- flexmap_to_long(parsed)

    ag_raw <- readr::read_csv(ag_path, show_col_types = FALSE)
    tr_raw <- readr::read_csv(tr_path, show_col_types = FALSE)

    va <- validate_antigen_table(ag_raw, long)
    vt <- validate_traceability_table(tr_raw, long)

    validate(
      need(va$ok, paste("Antigen table errors:", paste(va$errors, collapse = "; "))),
      need(vt$ok, paste("Traceability table errors:", paste(vt$errors, collapse = "; ")))
    )

    list(
      parsed = parsed, long = long,
      merged = merge_qc_data(long, va$data, vt$data),
      warnings = c(va$warnings, vt$warnings)
    )
  })

  output$data_status <- renderUI({
    d <- loaded()
    n_ok <- tags$p(class = "text-success",
                   sprintf("Loaded: %d wells, %d analytes.",
                          dplyr::n_distinct(d$merged$well_384), dplyr::n_distinct(d$merged$analyte)))
    if (length(d$warnings)) {
      tagList(n_ok, tags$p(class = "text-warning small", paste(d$warnings, collapse = " ")))
    } else n_ok
  })

  observe({
    choices <- plate_view_choices(loaded()$merged)
    updateSelectInput(session, "plate_view_mode", choices = choices, selected = "min")
  })

  observe({
    ants <- sort(unique(loaded()$merged$antigen))
    default <- if ("CHIKV_E1" %in% ants) "CHIKV_E1" else ants[1]
    updateSelectInput(session, "facet_antigen", choices = ants, selected = default)
  })

  thresholds_r <- reactive(mfi_control_thresholds(loaded()$merged))
  flags_r <- reactive(build_sample_flags(loaded()$merged, thresholds = thresholds_r()))

  # ---- Bead count QC ----
  output$plate_plot <- renderGirafe({
    plot_plate_beadcount(loaded()$merged, view = input$plate_view_mode)
  })

  # ---- MFI QC ----
  output$mfi_controls <- renderGirafe(plot_mfi_vs_controls(loaded()$merged))
  output$mfi_antigen   <- renderGirafe(plot_mfi_by_antigen(loaded()$merged))
  output$mfi_sample    <- renderGirafe(plot_mfi_by_sample(loaded()$merged, thresholds = thresholds_r()))
  output$mfi_cv        <- renderGirafe(plot_pool_cv(loaded()$merged))
  output$mfi_facet <- renderGirafe({
    req(input$facet_antigen)
    plot_mfi_by_antigen_facet(loaded()$merged, antigens = input$facet_antigen)
  })

  # ---- Refinement ----
  output$flags_table <- renderDT({
    flagged_for_review(flags_r()) |> datatable(options = list(pageLength = 10))
  })

  qc_result <- eventReactive(input$apply_removal, {
    manual_ids <- if (!is.null(input$manual_removal_csv)) {
      read_manual_removal_list(input$manual_removal_csv$datapath)
    } else character(0)
    apply_qc_removals(
      loaded()$merged, flags_r(),
      remove_bead_categories = if (is.null(input$remove_categories)) character(0) else input$remove_categories,
      remove_low_mfi = FALSE,
      manual_sample_ids = manual_ids
    )
  }, ignoreNULL = FALSE)

  output$removal_log <- renderText(paste(qc_result()$log, collapse = "\n"))

  clean_data <- reactive({
    res <- tryCatch(qc_result(), error = function(e) NULL)
    if (is.null(res)) loaded()$merged else res$clean
  })

  # ---- Export ----
  output$dl_antigen_template <- downloadHandler(
    filename = "antigen_table_template.csv",
    content = function(file) readr::write_csv(antigen_template(parse_flexmap_csv(DEFAULT_RAW)), file)
  )
  output$dl_trace_template <- downloadHandler(
    filename = "traceability_table_template.csv",
    content = function(file) {
      long0 <- flexmap_to_long(parse_flexmap_csv(DEFAULT_RAW))
      readr::write_csv(traceability_template(long0), file)
    }
  )
  output$dl_sample_removal_template <- downloadHandler(
    filename = "sample_removal_template.csv",
    content = function(file) readr::write_csv(sample_removal_template(loaded()$merged), file)
  )
  output$dl_clean <- downloadHandler(
    filename = "qc_cleaned_long.csv",
    content = function(file) export_clean_long(clean_data(), file)
  )
  output$dl_flags <- downloadHandler(
    filename = "qc_flags.csv",
    content = function(file) export_qc_flags(flagged_for_review(flags_r()), file)
  )
  output$dl_report <- downloadHandler(
    filename = "qc_report.pdf",
    contentType = "application/pdf",
    content = function(file) {
      res <- tryCatch(qc_result(), error = function(e) NULL)
      export_qc_report_pdf(
        file, loaded()$parsed, thresholds_r(), flags_r(), loaded()$merged,
        removal_log = if (!is.null(res)) res$log else NULL
      )
    }
  )
}

shinyApp(ui, server)

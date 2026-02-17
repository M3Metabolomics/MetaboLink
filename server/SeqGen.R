
###############################
# Hygge Project – Seq generator
###############################

# Reactive storage for uploaded / edited data
Data_in <- reactiveVal(NULL)

#========================================================
# 1. File upload
#========================================================

# Store uploaded CSV in Data_in
observeEvent(input$dataFile, {
  df <- read.csv(input$dataFile$datapath, stringsAsFactors = FALSE)
  Data_in(df)
})

#========================================================
# 2. Helpers – smart column guessing
#========================================================

guess_column <- function(cols,
                         current = NULL,
                         exact_candidates = character(),
                         pattern_candidates = character()) {
  # 1) Keep current selection if still valid
  if (!is.null(current) && current %in% cols) {
    return(current)
  }
  
  # 2) Try exact candidates (case-insensitive)
  if (length(exact_candidates)) {
    for (cand in exact_candidates) {
      idx <- which(tolower(cols) == tolower(cand))
      if (length(idx) > 0) {
        return(cols[idx[1]])
      }
    }
  }
  
  # 3) Try regex / partial matches
  if (length(pattern_candidates)) {
    for (pat in pattern_candidates) {
      hits <- grep(pat, cols, ignore.case = TRUE)
      if (length(hits) > 0) {
        return(cols[hits[1]])
      }
    }
  }
  
  # 4) Fallback: just use the first column
  cols[1]
}

#========================================================
# 3. Column mapping dropdowns (name / position / type)
#========================================================

observe({
  df <- Data_in()
  if (is.null(df)) return()
  
  cols <- names(df)
  
  # Keep current selections if possible
  current_name <- isolate(input$col_name_col)
  current_pos  <- isolate(input$col_position_col)
  current_type <- isolate(input$col_sampletype_col)
  
  # NAME column: typical name / sample name patterns
  guess_name <- guess_column(
    cols,
    current = current_name,
    exact_candidates = c("name", "Name", "sample_name", "Sample.Name", "Sample ID"),
    pattern_candidates = c("name", "sample.?id", "sample.?name")
  )
  
  # POSITION / VIAL column
  guess_pos <- guess_column(
    cols,
    current = current_pos,
    exact_candidates = c("Position", "Vial", "Well", "WellID"),
    pattern_candidates = c("position", "vial", "well", "plate.*pos")
  )
  
  # SAMPLE TYPE column
  guess_type <- guess_column(
    cols,
    current = current_type,
    exact_candidates = c("Sample.type", "SampleType", "Type", "Sample_Type"),
    pattern_candidates = c("sample.?type", "^type$", "class", "group")
  )
  
  updateSelectInput(
    session, "col_name_col",
    choices = cols,
    selected = guess_name
  )
  updateSelectInput(
    session, "col_position_col",
    choices = cols,
    selected = guess_pos
  )
  updateSelectInput(
    session, "col_sampletype_col",
    choices = cols,
    selected = guess_type
  )
  updateSelectInput(session, "col_remove_select", choices = cols)
})

#========================================================
# 4. Overview of sample types
#========================================================

output$hygge_overview <- renderUI({
  df <- Data_in()
  
  # No file loaded
  if (is.null(df)) {
    return(
      div(
        class = "hygge-subtle",
        "No file loaded yet. Upload a CSV to see an overview."
      )
    )
  }
  
  stype_col <- input$col_sampletype_col
  
  # Sample type column not set or invalid
  if (is.null(stype_col) || !stype_col %in% names(df)) {
    return(
      div(
        p(strong("Sample type column not set.")),
        p(
          class = "hygge-subtle",
          "Select the column containing Sample / Blank / QC in the dropdown under 'Step 1'."
        )
      )
    )
  }
  
  stype <- df[[stype_col]]
  
  # Counts
  n_total  <- nrow(df)
  n_blank  <- sum(tolower(stype) == "blank",  na.rm = TRUE)
  n_qc     <- sum(tolower(stype) == "qc",     na.rm = TRUE)
  n_sample <- sum(tolower(stype) == "sample", na.rm = TRUE)
  
  # Unknown / missing types
  recognized  <- tolower(stype) %in% c("blank", "qc", "sample")
  unknown_idx <- which(is.na(stype) | stype == "" | !recognized)
  n_unknown   <- length(unknown_idx)
  
  tagList(
    p(strong("Overview of loaded data")),
    tags$ul(
      tags$li(paste("Total rows:", n_total)),
      tags$li(paste("Samples:", n_sample)),
      tags$li(paste("Blanks:", n_blank)),
      tags$li(paste("QCs:", n_qc)),
      tags$li(paste("Rows with unknown / missing type:", n_unknown))
    ),
    if (n_unknown > 0) {
      div(
        p(strong("Rows with unknown type:")),
        p(
          paste(head(unknown_idx, 20), collapse = ", "),
          if (n_unknown > 20) " ..." else ""
        ),
        p(
          class = "hygge-subtle",
          "These are row indices in the uploaded table where the selected type column is not 'Sample', 'Blank' or 'QC'."
        )
      )
    } else {
      div(
        class = "hygge-subtle",
        "All rows have a recognized type (Sample, Blank or QC)."
      )
    }
  )
})

#========================================================
# 5. Uploaded data table + editing
#========================================================

output$uploadedTable <- DT::renderDT({
  req(Data_in())
  DT::datatable(
    Data_in(),
    editable  = TRUE,          # user can edit cells directly
    selection = "multiple",    # for row deletion
    options   = list(pageLength = 100),
    caption   = "Unprocessed data"
  )
})

# Apply cell edits from DT
observeEvent(input$uploadedTable_cell_edit, {
  df   <- Data_in()
  info <- input$uploadedTable_cell_edit
  i <- info$row
  j <- info$col
  v <- info$value
  
  df[i, j] <- DT::coerceValue(v, df[i, j])  # keep column type
  Data_in(df)
})

# Delete selected rows
observeEvent(input$delete_rows, {
  df <- Data_in()
  req(df)
  
  s <- input$uploadedTable_rows_selected
  if (length(s) > 0) {
    df <- df[-s, , drop = FALSE]
    Data_in(df)
  }
})

# Add a completely empty row at the bottom
observeEvent(input$add_empty_row, {
  df <- Data_in()
  req(df)
  
  # Create a one-row data.frame with the same columns and NA values
  new_row <- as.data.frame(
    lapply(df, function(col) {
      if (is.factor(col)) {
        factor(NA, levels = levels(col))
      } else {
        NA
      }
    }),
    stringsAsFactors = FALSE
  )
  
  df2 <- rbind(df, new_row)
  Data_in(df2)
})

# Duplicate selected rows and append them at the bottom
observeEvent(input$dup_rows_btn, {
  df <- Data_in()
  req(df)
  
  s <- input$uploadedTable_rows_selected
  if (length(s) == 0) {
    showNotification(
      "Select at least one row in 'Uploaded data' before duplicating.",
      type = "error"
    )
    return()
  }
  
  df2 <- rbind(df, df[s, , drop = FALSE])
  Data_in(df2)
})


#========================================================
# 6. Column add / remove tools
#========================================================

# Add new column and auto-fill it
observeEvent(input$add_new_col_btn, {
  df <- Data_in()
  req(df)
  
  new_name <- input$new_col_name
  new_val  <- input$new_col_val
  
  if (is.null(new_name) || new_name == "") {
    showNotification("Please enter a name for the new column.", type = "error")
    return()
  }
  
  df[[new_name]] <- new_val
  Data_in(df)
  
  updateTextInput(session, "new_col_name", value = "")
  updateTextInput(session, "new_col_val", value = "")
  showNotification(paste("Added column:", new_name), type = "message")
})

# Remove selected column
observeEvent(input$remove_col_btn, {
  df <- Data_in()
  req(df)
  
  col_to_remove <- input$col_remove_select
  req(col_to_remove)
  
  if (col_to_remove %in% names(df)) {
    df[[col_to_remove]] <- NULL
    Data_in(df)
    showNotification(paste("Removed column:", col_to_remove), type = "message")
  }
})



#========================================================
# 7. Sequence generation (data_processed)
#========================================================

data_processed <- eventReactive(input$process, {
  df <- Data_in()
  req(df)
  
  # Column mapping
  name_col  <- input$col_name_col
  pos_col   <- input$col_position_col
  stype_col <- input$col_sampletype_col
  req(name_col, pos_col, stype_col)
  
  if (!all(c(name_col, pos_col, stype_col) %in% names(df))) {
    stop("Selected columns for name/position/type are not found in the data.")
  }
  
  stype <- df[[stype_col]]
  
  # User parameters
  num_eqQCs            <- input$num_eqQCs
  num_MSMSs            <- input$num_MSMSs
  num_QCs              <- input$num_QCs
  insert_after_samples <- input$insert_after_samples
  
  # Split data into subsets (by selected sample-type column)
  blanks  <- df[tolower(stype) == "blank",  , drop = FALSE]
  qcs     <- df[tolower(stype) == "qc",     , drop = FALSE]
  samples <- df[tolower(stype) == "sample", , drop = FALSE]
  
  if (nrow(qcs) == 0) {
    stop("No QC samples available to duplicate.")
  }
  
  existing_qcs <- qcs
  
  # Base name templates for generated rows
  # eqQC and MSMS keep their specific tags
  sample_type_eqQC_names <- paste0(existing_qcs[[name_col]], "_eq_QC")
  sample_type_MSMS_names <- paste0(existing_qcs[[name_col]], "_MSMS")
  
  # --- FIX APPLIED HERE ---
  # Removed paste0(..., "_QC") so it doesn't double up (e.g. avoid RPposQC_QC05)
  # It now uses the name exactly as it appears in the QC row
  new_QC_names_base      <- existing_qcs[[name_col]]
  new_QC_samples_names   <- existing_qcs[[name_col]]
  # ------------------------
  
  #----- Generate eqQCs -----
  if (num_eqQCs > 0) {
    new_eqQCs <- qcs[rep(1, num_eqQCs), , drop = FALSE]
    new_eqQC_names <- paste(
      sample_type_eqQC_names,
      sprintf("%02d", seq_len(num_eqQCs)),
      sep = ""
    )
    new_eqQCs[[name_col]] <- new_eqQC_names
  } else {
    new_eqQCs <- qcs[0, , drop = FALSE]
  }
  
  #----- Generate MSMS -----
  if (num_MSMSs > 0) {
    new_MSMSs <- qcs[rep(1, num_MSMSs), , drop = FALSE]
    new_MSMS_names <- paste(
      sample_type_MSMS_names,
      sprintf("%02d", seq_len(num_MSMSs)),
      sep = ""
    )
    new_MSMSs[[name_col]] <- new_MSMS_names
  } else {
    new_MSMSs <- qcs[0, , drop = FALSE]
  }
  
  #----- Generate QCs before samples -----
  if (num_QCs > 0) {
    new_QCs <- qcs[rep(1, num_QCs), , drop = FALSE]
    new_QC_names <- paste(
      new_QC_names_base,
      sprintf("%02d", seq_len(num_QCs)),
      sep = ""
    )
    new_QCs[[name_col]] <- new_QC_names
  } else {
    new_QCs <- qcs[0, , drop = FALSE]
  }
  
  # Combine only the QC blocks that have rows
  pieces    <- list(new_eqQCs, new_MSMSs, new_QCs)
  non_empty <- pieces[sapply(pieces, nrow) > 0]
  if (length(non_empty) > 0) {
    qcs_combined <- do.call(rbind, non_empty)
  } else {
    qcs_combined <- qcs[0, , drop = FALSE]
  }
  
  #----- Randomize samples and intersperse QCs -----
  samples_randomized <- samples[sample(nrow(samples)), , drop = FALSE]
  
  # Skeleton with the same columns as df
  samples_with_qcs <- df[0, , drop = FALSE]
  qc_counter <- num_QCs + 1
  
  for (i in seq_len(nrow(samples_randomized))) {
    # Add one sample row
    samples_with_qcs <- rbind(
      samples_with_qcs,
      samples_randomized[i, , drop = FALSE]
    )
    
    # After every N samples, add a QC row
    if (insert_after_samples > 0 &&
        i %% insert_after_samples == 0) {
      
      qc_row <- qcs[1, , drop = FALSE]   # clone structure of a QC row
      qc_row[[name_col]]  <- paste(
        new_QC_samples_names,
        sprintf("%02d", qc_counter),
        sep = ""
      )
      qc_row[[stype_col]] <- "QC"
      
      samples_with_qcs <- rbind(samples_with_qcs, qc_row)
      qc_counter <- qc_counter + 1
    }
  }
  
  # Final combined dataset: blanks, generated QCs, and sample/QC block
  Data_out <- rbind(blanks, qcs_combined, samples_with_qcs)
  
  #----- Optional: force sequence to end with a QC -----
  if (isTRUE(input$end_with_qc)) {
    if (nrow(Data_out) > 0) {
      last_type <- Data_out[[stype_col]][nrow(Data_out)]
      if (is.na(last_type) || tolower(last_type) != "qc") {
        qc_row <- qcs[1, , drop = FALSE]
        qc_row[[name_col]]  <- paste(
          new_QC_samples_names,
          sprintf("%02d", qc_counter),
          sep = ""
        )
        qc_row[[stype_col]] <- "QC"
        Data_out <- rbind(Data_out, qc_row)
      }
    }
  }
  # --- FIX: Reset weird row numbers (25.1, 25.2, etc.) to clean 1, 2, 3... ---
  rownames(Data_out) <- NULL
  
  
  Data_out
})  
# When processing is done, switch to the "Processed sequence" tab
observeEvent(input$process, {
  updateTabsetPanel(session, "hygge_tabs", selected = "Processed sequence")
})

#========================================================
# 8. Processed sequence table
#========================================================

output$table <- DT::renderDT({
  req(data_processed())
  DT::datatable(
    data_processed(),
    options = list(pageLength = 100),
    caption = "Data after processing"
  )
})

#========================================================
# 9. Bruker-format export (matches your Excel template)
#========================================================

brunker_data <- reactive({
  req(data_processed())
  data <- data_processed()
  
  name_col <- input$col_name_col
  pos_col  <- input$col_position_col
  req(name_col, pos_col)
  
  if (!all(c(name_col, pos_col) %in% names(data))) {
    stop("Selected name/position columns are not present in processed data.")
  }
  
  # Default values taken from your Bruker_format.xlsx
  sep_default       <- "D:\\Methods\\Users\\JH\\pro methods\\JH\\LC methods\\vanquish_lipid11.5min.m?HyStar_LC"
  ms_default        <- "D:\\Methods\\Users\\JH\\HT\\final\\Lipidomics\\Lipidomics_pos_1pasef_msmsStep_30ev_k0range_large.m?OtofImpacTEMControl"
  data_path_default <- "D:\\Data\\JH\\LPF"
  
  out <- data.frame(
    Vial                = data[[pos_col]],
    `Sample ID`         = data[[name_col]],
    `Method Set`        = NA_character_,
    `Separation Method` = sep_default,
    `Injection Method`  = NA_character_,
    `MS Method`         = ms_default,
    `Volume [µl]`       = 0.1,
    `Data Path`         = data_path_default,
    check.names         = FALSE,
    stringsAsFactors    = FALSE
  )
  
  out
})

#========================================================
# 10. Reset / clear data
#========================================================

observeEvent(input$reset_all, {
  Data_in(NULL)
  
  # Reset dropdowns to empty
  updateSelectInput(session, "col_name_col",       choices = character(0))
  updateSelectInput(session, "col_position_col",   choices = character(0))
  updateSelectInput(session, "col_sampletype_col", choices = character(0))
  updateSelectInput(session, "col_remove_select",  choices = character(0))
  
  updateTabsetPanel(session, "hygge_tabs", selected = "Uploaded data")
  showNotification("Data reset.", type = "warning")
})

#========================================================
# 11. Download handlers
#========================================================

# Processed sequence as CSV
output$downloadData <- downloadHandler(
  filename = function() {
    paste("data-output-", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    req(data_processed())
    write.csv(data_processed(), file, row.names = FALSE)
  }
)

# Bruker-format as XLSX
output$downloadData_bruker <- downloadHandler(
  filename = function() {
    paste("data-output-", Sys.Date(), ".xlsx", sep = "")
  },
  content = function(file) {
    req(brunker_data())
    openxlsx::write.xlsx(brunker_data(), file, rowNames = FALSE)
  }
)

#========================================================
# 12. Open modal from main UI button
#========================================================

observeEvent(input$hygge_open, {
  showModal(hygge_modal())
})


#========================================================
# 14. Step 6 - Final Statistics & Runtime
#========================================================

output$final_stats_ui <- renderUI({
  # 1. Get the processed data safely
  df <- tryCatch(data_processed(), error = function(e) NULL)
  
  # 2. If data hasn't been processed yet, show a polite message
  if (is.null(df)) {
    return(
      div(
        class = "hygge-subtle",
        style = "margin-top: 10px;",
        "Please run 'Process data' (Step 5) to see the final statistics."
      )
    )
  }
  
  # 3. Get the column name used for Sample Type
  stype_col <- input$col_sampletype_col
  
  # Validation: Ensure column exists
  if (is.null(stype_col) || !stype_col %in% names(df)) {
    return(div(class = "text-danger", "Error: Sample type column missing in processed data."))
  }
  
  # 4. Calculate Counts for the PROCESSED data
  stypes   <- tolower(df[[stype_col]])
  n_total  <- nrow(df)
  n_sample <- sum(stypes == "sample", na.rm = TRUE)
  n_blank  <- sum(stypes == "blank",  na.rm = TRUE)
  # Note: Generated eqQCs and MSMSs are usually labeled "QC" in the type column 
  # based on the generator logic, so this captures all QCs.
  n_qc     <- sum(stypes == "qc",     na.rm = TRUE)
  
  # 5. Calculate Runtime
  t_min <- input$method_duration
  if (is.na(t_min) || t_min < 0) t_min <- 0
  
  total_minutes <- n_total * t_min
  
  # Time formatting
  days    <- floor(total_minutes / (24 * 60))
  rem_min <- total_minutes %% (24 * 60)
  hours   <- floor(rem_min / 60)
  minutes <- round(rem_min %% 60)
  
  time_str <- paste0(hours, "h ", minutes, "m")
  if (days > 0) {
    time_str <- paste0(days, "d ", time_str)
  }
  
  # 6. Render the List
  tagList(
    div(
      style = "margin-top: 15px;",
      
      # Section A: Sequence Composition
      p(strong("Generated Sequence Stats:")),
      tags$ul(
        tags$li(paste("Total rows:", n_total)),
        tags$li(paste("Samples:", n_sample)),
        tags$li(paste("Blanks:", n_blank)),
        tags$li(paste("Total QCs:", n_qc)),
        tags$small(class="hygge-subtle", "(Includes eqQC, MSMS, and standard QCs)")
      ),
      
      tags$hr(),
      
      # Section B: Runtime
      p(strong("Estimated Runtime:")),
      tags$ul(
        tags$li(
          span("Total Time: ", style="color: #555;"),
          strong(time_str, style = "color: #2c3e50; font-size: 1.1em;")
        )
      )
    )
  )
})

#========================================================
# 13. Modal UI definition
#========================================================

hygge_modal <- function() {
  modalDialog(
    title = "Hygge Project seq generator",
    size = "xl",
    easyClose = TRUE,
    footer = modalButton("Close"),
    div(
      class = "hygge-wrap",
      fluidRow(
        # LEFT: Sidebar with steps
        column(
          width = 4,
          class = "hygge-sidebar",
          
          # STEP 1: Upload & map
          div(
            class = "hygge-table-box",
            h4("Step 1 · Upload & map columns"),
            
            fileInput(
              "dataFile", "Choose CSV file",
              accept = c("text/csv", ".csv"),
              buttonLabel = "Browse...",
              placeholder = "No file selected",
              width = "100%"
            ),
            
            # Reset button
            div(
              style = "margin-bottom: 15px;",
              actionButton(
                "reset_all",
                "Clear / Reset Data",
                icon  = icon("refresh"),
                class = "btn-xs btn-warning"
              )
            ),
            
            selectInput("col_name_col",      "Column for sample name",  choices = NULL),
            selectInput("col_position_col",  "Column for position / vial", choices = NULL),
            selectInput("col_sampletype_col","Column for sample type",  choices = NULL),
            
            tags$small(
              class = "hygge-subtle",
              "Upload a CSV first, then choose which columns are name, position and type."
            )
          ),
          
          # STEP 2: Data overview
          div(
            class = "hygge-table-box",
            h4("Step 2 · Data overview"),
            uiOutput("hygge_overview")
          ),
          
          # STEP 3: Generator settings
          div(
            class = "hygge-table-box",
            h4("Step 3 · Generator settings"),
            fluidRow(
              column(
                6,
                numericInput("num_eqQCs", "Number of eqQCs", value = 1, min = 0)
              ),
              column(
                6,
                numericInput("num_MSMSs", "Number of MSMSs", value = 6, min = 0)
              )
            ),
            fluidRow(
              column(
                6,
                numericInput("num_QCs", "QCs before samples", value = 4, min = 0)
              ),
              column(
                6,
                numericInput("insert_after_samples",
                             "Insert QC after every N samples",
                             value = 5, min = 0)
              )
            ),
            checkboxInput(
              "end_with_qc",
              "Force sequence to end with a QC",
              value = TRUE
            )
          ),
          
          # STEP 4: Edit data
          # STEP 4: Edit data
          div(
            class = "hygge-table-box",
            
            # --- START CHANGE: Wrap content in details/summary ---
            tags$details(
              
              # The clickable Header
              tags$summary(
                h4(
                  "Step 4 · Edit table data ",
                  # This creates the smaller, subtle text
                  tags$span("(Click to show)", style = "font-size: 0.6em; color: #888; font-weight: normal; margin-left: 5px;"),
                  style = "display:inline; cursor:pointer"
                )
              ),
              
              # The Content to hide/show
              div(
                style = "margin-top: 15px;", # Add a little spacing when opened
                
                tags$small("Edits apply to the 'Uploaded data' tab."),
                tags$small("Psssst, double click in the data table, to edit any data."),
                
                br(), br(),
                
                # 4a. Add column
                strong("Add new column"),
                fluidRow(
                  column(6, textInput("new_col_name", NULL, placeholder = "Column name")),
                  column(6, textInput("new_col_val",  NULL, placeholder = "Fill value"))
                ),
                actionButton(
                  "add_new_col_btn",
                  "Add & Fill",
                  icon  = icon("plus"),
                  class = "btn-xs"
                ),
                br(), br(),
                
                # 4b. Remove column
                strong("Remove column"),
                fluidRow(
                  column(8, selectInput("col_remove_select", NULL, choices = NULL)),
                  column(
                    4,
                    actionButton(
                      "remove_col_btn",
                      "Remove",
                      icon  = icon("trash"),
                      class = "btn-xs btn-danger"
                    )
                  )
                ),
                br(), br(),
                
                # 4c. Delete rows
                strong("Delete row(s)"),
                tags$small(
                  class = "hygge-subtle",
                  "Select one or more rows in 'Uploaded data', then click delete."
                ),
                br(),
                actionButton(
                  "delete_rows",
                  "Delete selected rows from table",
                  class = "btn-xs"
                ),
                
                br(), br(),
                
                # 4d. Add / duplicate rows
                strong("Add row(s)"),
                tags$small(
                  class = "hygge-subtle",
                  "Append an empty row or duplicate selected rows."
                ),
                br(),
                actionButton(
                  "add_empty_row",
                  "Add empty row",
                  class = "btn-xs"
                ),
                actionButton(
                  "dup_rows_btn",
                  "Duplicate selected rows",
                  class = "btn-xs"
                )
              ) 
            ) 
            # --- END CHANGE ---
          ),
          
          # STEP 5: Run & download
          div(
            class = "hygge-table-box",
            h4("Step 5 · Run & download"),
            
            strong("1. Process sequence"),
            br(),
            bsButton("process", "Process data", style = "primary"),
            br(), br(),
            
            strong("2. Download output"),
            div(
              class = "hygge-inline-buttons",
              downloadButton("downloadData_bruker", "Bruker (.xlsx)"),
              downloadButton("downloadData",        "Generic (.csv)")
            )
          ),
          
          # STEP 6: Final Overview & Runtime
          div(
            class = "hygge-table-box",
            h4("Step 6 · Final Overview & Runtime"),
            
            # Input for method duration
            numericInput(
              "method_duration", 
              "Method duration inc. injection (min)", 
              value = 15, 
              min = 0.1, 
              step = 0.5
            ),
            
            # The Combined Output (Stats + Runtime)
            uiOutput("final_stats_ui")
          )
          
        ),
        
        # RIGHT: Tabs – uploaded vs processed
        column(
          width = 8,
          div(
            class = "hygge-table-box hygge-dt",
            tabsetPanel(
              id = "hygge_tabs",
              tabPanel(
                title = "Uploaded data",
                h5("Uploaded CSV (source)"),
                DT::DTOutput("uploadedTable")
              ),
              tabPanel(
                title = "Processed sequence",
                h5("Generated sequence (result)"),
                DT::DTOutput("table")
              )
            )
          )
        )
      )
    )
  )
}




# ==============================================================================
# SERVER: AQ logic
# ==============================================================================

# ------------------------------------------------------------------------------
# UI builder for Absolute Quantification modal
# ------------------------------------------------------------------------------
absq_modal_body <- function() {
  tagList(
    div(
      class = "absq-modal-body",
      
      # 1. Settings -------------------------------------------------------------
      div(
        class = "absq-section",
        h4("1. Settings"),
        fluidRow(
          column(
            6,
            radioButtons(
              "absq_matrix_type", "Sample Matrix",
              choices  = c("Tissue (mg)", "Biofluid (µL)"),
              selected = "Tissue (mg)",
              width    = "100%"
            )
          ),
          column(
            6,
            checkboxInput(
              "absq_include_qc",
              "Include QC samples in checks",
              value = TRUE,
              width = "100%"
            )
          )
        )
      ),
      
      # 2. Readiness checks ----------------------------------------------------
      div(
        class = "absq-section",
        h4("2. Readiness Checks"),
        uiOutput("absq_readiness_ui")
      ),
      
      # 3. Internal standard spikes -------------------------------------------
      div(
        class = "absq-section",
        h4("3. Internal Standard Spikes"),
        
        fluidRow(
          column(
            6,
            selectInput(
              "absq_spike_unit",
              "Input Unit",
              choices  = c("pmol", "nmol"),
              selected = "pmol",
              width    = "100%"
            )
          ),
          column(
            6,
            div(
              class = "absq-help-text",
              "Double-click cells in the table to edit individual amounts."
            )
          )
        ),
        
        br(),
        
        # Auto-fill panel for spikes
        wellPanel(
          style = "padding: 10px; background-color: #fff; border: 1px solid #ddd;",
          p(
            style = "font-size: 0.9em; color: #666; margin-bottom: 5px;",
            tags$i(class = "fa fa-info-circle"),
            " Psst, if you have the same value for all of 'spiked amount', just auto fill here:"
          ),
          fluidRow(
            column(
              6,
              numericInput(
                "absq_autofill_val",
                "Auto-fill Amount:",
                value = 0,
                step  = 0.1,
                width = "100%"
              )
            ),
            column(
              6,
              style = "margin-top: 25px;",
              actionButton(
                "absq_autofill_btn",
                "Apply to All Rows",
                icon  = icon("arrow-down"),
                class = "btn-warning",
                width = "100%"
              )
            )
          )
        ),
        
        DT::DTOutput("absq_spike_table")
      ),
      
      # 4. Sample Amounts (Updated UI) ----------------------------------
      div(
        class = "absq-section",
        h4("4. Sample Amounts"),
        p(
          "Provide the sample mass / volume for each sample. ",
          "If available, values are pre-filled from the sequence file, ",
          "but you can override or enter them manually."
        ),
        
        wellPanel(
          style = "padding: 10px; background-color: #fff; border: 1px solid #ddd;",
          p(
            style = "font-size: 0.9em; color: #666; margin-bottom: 5px;",
            tags$i(class = "fa fa-info-circle"),
            " Bulk Actions:"
          ),
          fluidRow(
            # --- Existing Auto-fill Input ---
            column(
              4,
              numericInput(
                "absq_amt_autofill_val",
                "Auto-fill Value:",
                value = NA,
                step  = 0.1,
                width = "100%"
              )
            ),
            # --- Existing Apply All Button ---
            column(
              4,
              style = "margin-top: 25px;",
              actionButton(
                "absq_amt_autofill_btn",
                "Apply to All",
                icon  = icon("arrow-down"),
                class = "btn-info",
                width = "100%"
              )
            ),
            # --- NEW BUTTON: Impute Average ---
            column(
              4, 
              style = "margin-top: 25px;",
              tipify(
                actionButton(
                  "absq_impute_avg_btn",
                  "Fill NAs with Avg",
                  icon  = icon("magic"),
                  class = "btn-success", # Green button to stand out
                  width = "100%"
                ),
                title = "Calculates average of existing values and fills empty cells only."
              )
            )
          )
        ),
        
        DT::DTOutput("absq_amount_table")
      ),
      
      # 5. Compute & results ---------------------------------------------------
      div(
        class = "absq-section",
        h4("5. Compute & Results"),
        
        selectInput(
          "absq_is_method",
          "IS Matching Strategy",
          choices  = c("Nearest RT", "Same lipid structure"),
          selected = "Nearest RT",
          width    = "50%"
        ),
        
        bsButton(
          "absq_compute",
          "Compute Absolute Quantification",
          style = "primary",
          width = "100%"
        ),
        
        tags$hr(),
        h5("Preview Results"),
        DT::DTOutput("absq_result_table_modal"),
        
        
        tags$hr(),
        checkboxInput(
          "absq_save_as_new",
          "Save as new dataset",
          value = TRUE,
          width = "100%"
        )
      )
    )
  )
}

# ==============================================================================
# Absolute Quantification (AQ) – SERVER BLOCK
# ==============================================================================

# ---------------------------------
# Helper functions
# ---------------------------------

# Bullet with green check / red cross for readiness UI
absq_readiness_bullet <- function(ok, txt) {
  color <- if (ok) "#2e7d32" else "#c62828"
  icon  <- shiny::icon(if (ok) "check" else "times")
  
  shiny::tags$li(
    style = paste0("color:", color, "; font-weight:600;"),
    icon,
    shiny::span(style = "color:black; font-weight:normal; margin-left:4px;", txt)
  )
}

# Find an 'amount' column in sequence table (exact match preferred, then fuzzy)
absq_find_amount_col <- function(seq_tbl) {
  exact <- which(tolower(colnames(seq_tbl)) == "amount")
  if (length(exact) >= 1) return(exact[1])
  
  fuzzy <- grep("amount", colnames(seq_tbl), ignore.case = TRUE)
  if (length(fuzzy) >= 1) return(fuzzy[1])
  
  NA_integer_
}

# Find RT column index, either via sequence labels or data column names
absq_find_rt_col <- function(seq_tbl, d) {
  rt_from_seq <- which(seq_tbl[, 1] == "RT")
  if (length(rt_from_seq)) return(rt_from_seq[1])
  
  rt_from_data <- which(colnames(d) == "RT")
  if (length(rt_from_data)) return(rt_from_data[1])
  
  NA_integer_
}

# Choose IS index for each feature (Nearest RT / Same lipid structure)
absq_choose_is_index <- function(d, seq_tbl, is_rows,
                                 method = c("Nearest RT", "Same lipid structure")) {
  method <- match.arg(method)
  
  rt_col_idx   <- absq_find_rt_col(seq_tbl, d)
  name_rows    <- which(seq_tbl[, 1] == "Name")
  name_col_idx <- if (length(name_rows)) name_rows[1] else 1
  
  n_feat <- nrow(d)
  
  # Base mapping: nearest RT
  if (!is.na(rt_col_idx) && length(is_rows) > 0) {
    feat_rts <- suppressWarnings(as.numeric(d[, rt_col_idx]))
    is_rts   <- suppressWarnings(as.numeric(d[is_rows, rt_col_idx]))
    
    nearest <- vapply(
      feat_rts,
      function(y) {
        if (is.na(y)) return(1L)
        which.min(abs(is_rts - y))
      },
      integer(1)
    )
  } else {
    nearest <- rep(1L, n_feat) # fallback
  }
  
  if (method == "Nearest RT") return(nearest)
  
  # Same lipid structure: match by lipid class prefix
  feat_names <- as.character(d[, name_col_idx])
  is_names   <- as.character(d[is_rows, name_col_idx])
  
  feat_class <- sub("^([A-Za-z]+).*", "\\1", feat_names)
  is_class   <- sub("^([A-Za-z]+).*", "\\1", is_names)
  
  chosen <- nearest
  
  if (!is.na(rt_col_idx) && length(is_rows) > 0) {
    is_rts   <- suppressWarnings(as.numeric(d[is_rows, rt_col_idx]))
    feat_rts <- suppressWarnings(as.numeric(d[, rt_col_idx]))
  }
  
  for (i in seq_len(n_feat)) {
    matches <- which(is_class == feat_class[i])
    if (length(matches) > 0) {
      if (!is.na(rt_col_idx)) {
        current_rt <- feat_rts[i]
        if (!is.na(current_rt)) {
          local_rts <- is_rts[matches]
          best_idx  <- which.min(abs(local_rts - current_rt))
          chosen[i] <- matches[best_idx]
        } else {
          chosen[i] <- matches[1]
        }
      } else {
        chosen[i] <- matches[1]
      }
    }
  }
  chosen
}

# Sticky footer for modal buttons
sticky_footer <- function(...) {
  shiny::tagList(
    tags$div(
      ...,
      style = paste0(
        "position: fixed;",
        "bottom: 0;",
        "left: 0;",
        "width: 100%;",
        "padding: 10px 15px;",
        "border-top: 1px solid #e5e5e5;",
        "background-color: #f5f5f5;",
        "z-index: 1050;"
      )
    )
  )
}

# ------------------------------------------------------------------------------
# Basic data access reactives
# ------------------------------------------------------------------------------

absq_data <- reactive({
  req(rv$activeFile)
  rv$data[[rv$activeFile]]
})

absq_seq <- reactive({
  req(rv$activeFile)
  rv$sequence[[rv$activeFile]]
})

# Which sample columns are “in play”
absq_selected_sample_mask <- reactive({
  req(absq_seq(), !is.null(input$absq_include_qc))
  sel_labels <- if (isTRUE(input$absq_include_qc)) c("Sample", "QC") else "Sample"
  absq_seq()[, 1] %in% sel_labels
})

# Internal standards (rows + names)
absq_is_parsed <- reactive({
  d <- absq_data()
  is_idx <- grep("\\(IS\\)", d[, 1], ignore.case = TRUE)
  if (length(is_idx) == 0) return(NULL)
  
  data.frame(
    IS_row  = is_idx,
    IS_name = d[is_idx, 1],
    stringsAsFactors = FALSE
  )
})

# ------------------------------------------------------------------------------
# State objects
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Auto-Initialization Logic (FIXED)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Auto-Initialization Logic (FIXED)
# ------------------------------------------------------------------------------

# This observer ensures rv$absq_amounts is ALWAYS populated correctly 
# whenever the modal is open or sample selection changes.
observe({
  req(absq_data(), absq_seq(), absq_selected_sample_mask())
  
  d   <- absq_data()
  seq <- absq_seq()
  
  mask        <- absq_selected_sample_mask()
  sample_cols <- which(mask)
  if (!length(sample_cols)) return()
  
  sample_ids <- colnames(d)[sample_cols]
  
  # 1. Check if the Sequence File actually has valid data
  seq_has_data <- FALSE
  amt_col_idx  <- absq_find_amount_col(seq)
  if (!is.na(amt_col_idx)) {
    # Check if there are non-NA values for the selected samples
    seq_vals <- suppressWarnings(as.numeric(seq[sample_ids, amt_col_idx]))
    if (any(!is.na(seq_vals))) seq_has_data <- TRUE
  }
  
  # 2. Check current internal table state
  # === THE FIX IS BELOW ===
  # We use isolate() so this observer does NOT trigger when the user edits the table
  curr <- isolate(rv$absq_amounts) 
  # ========================
  
  internal_has_nas <- FALSE
  if (!is.null(curr)) {
    if (any(is.na(curr$Amount))) internal_has_nas <- TRUE
  }
  
  # 3. Decision: Should we rebuild the internal table?
  rebuild <- FALSE
  if (is.null(curr)) {
    # Case A: Table doesn't exist yet -> Rebuild
    rebuild <- TRUE
  } else if (!identical(curr$Sample, sample_ids)) {
    # Case B: Samples selected have changed -> Rebuild
    rebuild <- TRUE
  } else if (internal_has_nas && seq_has_data) {
    # Case C: Internal table has NAs, but Sequence file HAS data.
    # This forces a sync to pull the values from the sequence file.
    rebuild <- TRUE
  }
  
  if (rebuild) {
    amt <- rep(NA_real_, length(sample_ids))
    if (!is.na(amt_col_idx)) {
      vals <- suppressWarnings(as.numeric(seq[sample_ids, amt_col_idx]))
      amt  <- vals
    }
    
    rv$absq_amounts <- data.frame(
      Sample = sample_ids,
      Amount = amt,
      stringsAsFactors = FALSE
    )
  }
})

# Check if sequence file has valid amounts (Status check)
absq_amounts_ready <- reactive({
  req(absq_seq())
  seq <- absq_seq()
  has_amt_col <- "amount" %in% tolower(colnames(seq))
  if (!has_amt_col) return(FALSE)
  
  # Are there any non-NA values?
  amt_col <- absq_find_amount_col(seq)
  if (is.na(amt_col)) return(FALSE)
  
  # Check if at least one selected sample has a value
  mask <- tryCatch(absq_selected_sample_mask(), error = function(e) TRUE)
  if (is.logical(mask) && length(mask) == nrow(seq)) {
    # Use only selected samples for validity check
    vals <- seq[mask, amt_col]
  } else {
    vals <- seq[, amt_col]
  }
  
  any(!is.na(suppressWarnings(as.numeric(vals))))
})

# ------------------------------------------------------------------------------
# Open modal
# ------------------------------------------------------------------------------

observeEvent(input$absq_open_modal, {
  # Open modal immediately
  showModal(
    modalDialog(
      title  = "Absolute Quantification",
      size   = "l",
      footer = sticky_footer(
        modalButton("Close window"),
        bsButton("absq_save", "Save & Apply",
                 style = "success", icon = icon("save"))
      ),
      absq_modal_body()
    )
  )
  
  # Reset result state
  rv$absq_result <- NULL
  shinyjs::disable("absq_save")
})

# ------------------------------------------------------------------------------
# Readiness checks (depends on matrix type)
# ------------------------------------------------------------------------------

output$absq_readiness_ui <- renderUI({
  req(input$absq_matrix_type, absq_seq(), absq_data())
  
  seq     <- absq_seq()
  d       <- absq_data()
  is_info <- absq_is_parsed()
  
  has_amt_col <- "amount" %in% tolower(colnames(seq))
  
  # Check 1: Sequence has amounts?
  amt_seq_ok <- absq_amounts_ready()
  
  # Check 2: Manual input has amounts?
  amt_manual_ok <- !is.null(rv$absq_amounts) &&
    any(!is.na(suppressWarnings(as.numeric(rv$absq_amounts$Amount))))
  
  # Combined check
  amt_any_ok <- amt_seq_ok || amt_manual_ok
  
  has_is  <- !is.null(is_info) && nrow(is_info) > 0
  n_is    <- if (is.null(is_info)) 0 else nrow(is_info)
  
  has_rt  <- "RT" %in% seq[, 1] || "RT" %in% colnames(d)
  
  is_bio   <- grepl("Biofluid", input$absq_matrix_type, ignore.case = TRUE)
  unit_txt <- if (is_bio) "volume in µL" else "mass in mg"
  
  tags$ul(
    absq_readiness_bullet(
      has_amt_col,
      paste0("Sequence contains 'amount' column (", unit_txt, ")")
    ),
    absq_readiness_bullet(
      amt_any_ok,
      "Sample amounts available (from sequence or manual input)"
    ),
    absq_readiness_bullet(
      has_is,
      paste0("Internal Standards detected (", n_is, " found)")
    ),
    absq_readiness_bullet(
      has_rt,
      "Retention Time (RT) available"
    )
  )
})

# ------------------------------------------------------------------------------
# Dynamic UI: Manual Sample Amounts
# ------------------------------------------------------------------------------

output$absq_manual_amounts_ui <- renderUI({
  # Only show this section if readiness check for amounts is NEGATIVE
  # (i.e. sequence doesn't have amounts)
  ready <- absq_amounts_ready()
  if (isTRUE(ready)) return(NULL) # Hide if we have amounts
  
  div(
    class = "absq-section",
    h4("Missing Sample Amounts"),
    p(
      "The sequence file does not contain sample amounts. ",
      "Please enter the sample mass / volume manually below."
    ),
    
    wellPanel(
      style = "padding: 10px; background-color: #fff; border: 1px solid #ddd;",
      p(
        style = "font-size: 0.9em; color: #666; margin-bottom: 5px;",
        tags$i(class = "fa fa-info-circle"),
        " If all samples use the same amount, you can auto-fill that value here:"
      ),
      fluidRow(
        column(
          6,
          numericInput(
            "absq_amt_autofill_val",
            "Auto-fill Spiked Amount:",
            value = NA,
            step  = 0.1,
            width = "100%"
          )
        ),
        column(
          6,
          style = "margin-top: 25px;",
          actionButton(
            "absq_amt_autofill_btn",
            "Apply to All Samples",
            icon  = icon("arrow-down"),
            class = "btn-info",
            width = "100%"
          )
        )
      )
    ),
    
    DT::DTOutput("absq_amount_table")
  )
})


# ------------------------------------------------------------------
# NEW: Impute missing (NA) amounts with the average of existing ones
# ------------------------------------------------------------------

observeEvent(input$absq_impute_avg_btn, {
  req(rv$absq_amounts)
  
  # 1. Get current values
  df <- rv$absq_amounts
  vals <- suppressWarnings(as.numeric(df$Amount))
  
  # 2. Check if we have enough data to calculate an average
  valid_vals <- vals[!is.na(vals)]
  
  if (length(valid_vals) == 0) {
    showNotification("Cannot calculate average: No valid amounts entered yet.", type = "warning")
    return()
  }
  
  # 3. Identify missing indices
  missing_idx <- which(is.na(vals))
  
  if (length(missing_idx) == 0) {
    showNotification("No missing values (NA) found to fill.", type = "message")
    return()
  }
  
  # 4. Calculate Average and rounded
  avg_val <- mean(valid_vals)
  avg_display <- round(avg_val, 2)
  
  # 5. Update only the missing rows
  df$Amount[missing_idx] <- avg_display
  rv$absq_amounts <- df
  
  # 6. Success message
  msg <- paste0("Filled ", length(missing_idx), " missing samples with average (", avg_display, ")")
  
  # FIXED: Changed type="success" to type="message" to avoid crash
  showNotification(msg, type = "message", duration = 4) 
})

# ------------------------------------------------------------------------------
# Spike tables per matrix type
# ------------------------------------------------------------------------------

absq_ensure_spike_table <- function() {
  info <- absq_is_parsed()
  key  <- isolate(input$absq_matrix_type)
  if (is.null(info) || is.null(key)) return()
  
  if (is.null(rv$absq_spikes[[key]])) {
    base_df <- data.frame(
      IS_row     = info$IS_row,
      IS_name    = info$IS_name,
      spike_pmol = 0,
      stringsAsFactors = FALSE
    )
    
    # copy from other matrix type if available
    other_keys <- setdiff(c("Tissue (mg)", "Biofluid (µL)"), key)
    donor <- NULL
    for (ok in other_keys) {
      if (!is.null(rv$absq_spikes[[ok]])) {
        donor <- rv$absq_spikes[[ok]]
        break
      }
    }
    if (!is.null(donor)) {
      merged <- merge(
        base_df,
        donor[, c("IS_row", "spike_pmol")],
        by = "IS_row", all.x = TRUE, suffixes = c("", ".donor")
      )
      base_df$spike_pmol <- ifelse(
        is.na(merged$spike_pmol.donor),
        base_df$spike_pmol,
        merged$spike_pmol.donor
      )
    }
    
    base_df <- base_df[order(base_df$IS_row), ]
    rv$absq_spikes[[key]] <- base_df
  }
}

observeEvent(list(absq_is_parsed(), input$absq_matrix_type), {
  req(absq_is_parsed(), input$absq_matrix_type)
  absq_ensure_spike_table()
})

absq_spike_display <- reactive({
  req(absq_seq(), absq_data(), absq_is_parsed(), input$absq_matrix_type)
  
  key <- input$absq_matrix_type
  absq_ensure_spike_table()
  df <- rv$absq_spikes[[key]]
  req(df)
  
  seq_tbl <- absq_seq()
  d       <- absq_data()
  
  rt_col_idx <- which(seq_tbl[, 1] == "RT")[1]
  if (is.na(rt_col_idx)) rt_col_idx <- which(colnames(d) == "RT")[1]
  
  rt_vals <- if (!is.na(rt_col_idx)) {
    round(as.numeric(d[df$IS_row, rt_col_idx]), 2)
  } else {
    rep(NA_real_, nrow(df))
  }
  
  unit_scale <- if (identical(input$absq_spike_unit, "nmol")) 1/1000 else 1
  col_name   <- paste0("Spiked Amount (", input$absq_spike_unit, ")")
  
  out <- data.frame(
    `IS row`  = df$IS_row,
    `IS name` = df$IS_name,
    `RT`      = rt_vals,
    check.names = FALSE
  )
  out[[col_name]] <- round(df$spike_pmol * unit_scale, 6)
  out
})

# When matrix type or spike unit changes, clear results + disable save
observeEvent(list(input$absq_matrix_type, input$absq_spike_unit), {
  rv$absq_result <- NULL
  shinyjs::disable("absq_save")
})

# Auto-switch spike input unit based on matrix type
observeEvent(input$absq_matrix_type, {
  req(input$absq_matrix_type)
  if (grepl("Biofluid", input$absq_matrix_type)) {
    updateSelectInput(session, "absq_spike_unit", selected = "pmol")
  } else {
    updateSelectInput(session, "absq_spike_unit", selected = "nmol")
  }
})

# Auto-fill spike amounts
observeEvent(input$absq_autofill_btn, {
  req(input$absq_matrix_type)
  key <- input$absq_matrix_type
  absq_ensure_spike_table()
  req(rv$absq_spikes[[key]])
  
  val <- input$absq_autofill_val
  if (is.na(val)) {
    showNotification("Please enter a valid numeric amount.", type = "warning")
    return()
  }
  
  unit_scale <- if (identical(input$absq_spike_unit, "nmol")) 1000 else 1
  fill_pmol  <- val * unit_scale
  
  rv$absq_spikes[[key]]$spike_pmol <- fill_pmol
  showNotification("All spike amounts updated for this matrix type.", type = "message")
})

output$absq_spike_table <- DT::renderDT({
  req(absq_spike_display())
  DT::datatable(
    absq_spike_display(),
    editable = list(
      target  = "cell",
      disable = list(columns = 0:2)
    ),
    options   = list(pageLength = 10, dom = "t", scrollX = FALSE),
    selection = "none",
    rownames  = FALSE
  )
})

observeEvent(input$absq_spike_table_cell_edit, {
  info <- input$absq_spike_table_cell_edit
  req(info, input$absq_matrix_type)
  
  key <- input$absq_matrix_type
  absq_ensure_spike_table()
  req(rv$absq_spikes[[key]])
  
  val <- suppressWarnings(as.numeric(info$value))
  if (is.na(val)) return()
  
  unit_scale <- if (identical(input$absq_spike_unit, "nmol")) 1000 else 1
  
  disp    <- absq_spike_display()
  is_row  <- disp[info$row, "IS row"]
  idx_spk <- match(is_row, rv$absq_spikes[[key]]$IS_row)
  if (is.na(idx_spk)) return()
  
  rv$absq_spikes[[key]]$spike_pmol[idx_spk] <- val * unit_scale
})

# ------------------------------------------------------------------------------
# Sample amounts table (manual / fallback)
# ------------------------------------------------------------------------------

absq_amount_display <- reactive({
  # NOTE: We use req(rv$absq_amounts) here.
  # Since we added the observer to auto-init rv$absq_amounts, this should now be safe.
  req(rv$absq_amounts, input$absq_matrix_type)
  
  df <- rv$absq_amounts
  
  is_bio   <- grepl("Biofluid", input$absq_matrix_type, ignore.case = TRUE)
  unit_den <- if (is_bio) "\u00B5L" else "mg"
  
  df$Unit <- unit_den
  df
})

output$absq_amount_table <- DT::renderDT({
  req(absq_amount_display())
  
  DT::datatable(
    absq_amount_display(),
    editable = list(
      target  = "cell",
      disable = list(columns = c(0, 2))  # Sample + Unit read-only
    ),
    selection = "none",
    rownames  = FALSE,
    options   = list(
      dom            = "t",      # Only show table (no search box/info text)
      scrollY        = "400px",  # <--- Fixes height & enables vertical scrolling
      scrollCollapse = TRUE,     # <--- Nice formatting if list is short
      paging         = FALSE     # <--- Disables pagination (shows ALL rows)
    )
  )
})

observeEvent(input$absq_amt_autofill_btn, {
  req(rv$absq_amounts)
  
  val <- input$absq_amt_autofill_val
  if (is.null(val) || is.na(val)) {
    showNotification("Please enter a valid amount to auto-fill.", type = "warning")
    return()
  }
  
  rv$absq_amounts$Amount <- as.numeric(val)
})

observeEvent(input$absq_amount_table_cell_edit, {
  info <- input$absq_amount_table_cell_edit
  req(info, rv$absq_amounts)
  
  row   <- info$row
  value <- suppressWarnings(as.numeric(info$value))
  if (is.na(value)) return()
  
  rv$absq_amounts$Amount[row] <- value
})

# ------------------------------------------------------------------------------
# Compute absolute amounts
# ------------------------------------------------------------------------------

observeEvent(input$absq_compute, {
  req(input$absq_matrix_type, absq_seq(), absq_data(), absq_is_parsed())
  
  key <- input$absq_matrix_type
  absq_ensure_spike_table()
  spikes <- rv$absq_spikes[[key]]
  req(spikes)
  
  d   <- absq_data()
  seq <- absq_seq()
  
  # 1) sample columns
  mask        <- absq_selected_sample_mask()
  sample_cols <- which(mask)
  req(length(sample_cols) > 0)
  sample_ids  <- colnames(d)[sample_cols]
  
  # 2) amounts from manual/auto-filled table
  # NOTE: Even if the manual table is hidden, rv$absq_amounts exists and holds the values.
  
  if (is.null(rv$absq_amounts)) {
    showNotification(
      "Internal Error: Sample amounts table is not initialized.",
      type = "error"
    )
    return()
  }
  
  amt_tbl <- rv$absq_amounts
  
  # Check synchronization
  if (!all(sample_ids %in% amt_tbl$Sample)) {
    showNotification(
      "Sample amounts are out of sync. Please close and reopen the window.",
      type = "error"
    )
    return()
  }
  
  idx     <- match(sample_ids, amt_tbl$Sample)
  raw_amt <- suppressWarnings(as.numeric(amt_tbl$Amount[idx]))
  
  # Validation: Do we have valid amounts?
  if (any(is.na(raw_amt))) {
    showNotification(
      "Some sample amounts are missing (NA). Please check the 'Readiness' or 'Manual Amounts' section.",
      type = "error"
    )
    return()
  }
  
  is_bio   <- grepl("Biofluid", key, ignore.case = TRUE)
  unit_den <- if (is_bio) "\u00B5L" else "mg"
  unit_lbl <- paste0("pmol/", unit_den)
  
  rv$absq_amounts_used <- data.frame(
    Sample = sample_ids,
    Amount = raw_amt,
    Unit   = unit_den,
    stringsAsFactors = FALSE
  )
  
  # 3) match IS rows
  is_rows <- spikes$IS_row
  is_map  <- absq_choose_is_index(d, seq, is_rows, method = input$absq_is_method)
  
  # 4) math
  raw_int <- d[, sample_cols, drop = FALSE]
  is_int  <- d[is_rows[is_map], sample_cols, drop = FALSE]
  
  spike_amt_pmol <- spikes$spike_pmol[is_map]  # base = pmol
  
  eps       <- 1e-12
  res       <- (raw_int / pmax(is_int, eps)) * spike_amt_pmol
  res_final <- sweep(res, 2, raw_amt, "/")  # divide by mass/volume
  
  # 5) preview table
  name_row_idx <- which(seq[, 1] == "Name")[1]
  if (is.na(name_row_idx)) name_row_idx <- 1
  
  compounds  <- d[, name_row_idx]
  is_used_nm <- d[is_rows[is_map], name_row_idx]
  
  preview_df <- data.frame(
    Compound  = compounds,
    `IS used` = is_used_nm,
    Units     = unit_lbl,
    check.names = FALSE
  )
  preview_df <- cbind(preview_df, as.data.frame(res_final))
  
  rv$absq_result <- preview_df
  showNotification(paste("Absolute amounts calculated in", unit_lbl), type = "message")
  shinyjs::enable("absq_save")
})

# ------------------------------------------------------------------------------
# Render result tables
# ------------------------------------------------------------------------------

output$absq_result_table_modal <- DT::renderDT({
  req(rv$absq_result)
  DT::datatable(
    rv$absq_result,
    options = list(pageLength = 10, scrollX = TRUE),
    rownames = FALSE
  )
})

output$absq_amounts_table_modal <- DT::renderDT({
  req(rv$absq_amounts_used)
  DT::datatable(
    rv$absq_amounts_used,
    options = list(pageLength = 10, dom = "t"),
    rownames = FALSE
  )
})

# ------------------------------------------------------------------------------
# Save back to dataset
# ------------------------------------------------------------------------------

observeEvent(input$absq_save, {
  req(rv$absq_result, absq_data(), absq_seq(), input$absq_matrix_type)
  
  d   <- absq_data()
  res <- rv$absq_result
  
  meta_cols       <- c("Compound", "IS used", "Units")
  sample_cols_res <- setdiff(colnames(res), meta_cols)
  
  common <- intersect(colnames(d), sample_cols_res)
  if (!length(common)) {
    showNotification(
      "No matching sample columns between data and AQ result.",
      type = "error"
    )
    return()
  }
  
  d_out  <- d
  idx_d  <- match(common, colnames(d))
  idx_res <- match(common, colnames(res))
  d_out[, idx_d] <- as.matrix(res[, idx_res, drop = FALSE])
  
  is_bio   <- grepl("Biofluid", input$absq_matrix_type, ignore.case = TRUE)
  unit_den <- if (is_bio) "uL" else "mg"
  suffix   <- paste0("_AbsQ_", unit_den)
  info     <- paste0("Absolute Quant (", input$absq_is_method, ", ", input$absq_spike_unit, ")")
  
  rv$tmpData     <- d_out
  rv$tmpSequence <- absq_seq()
  
  updateDataAndSequence(
    "Absolute Quantification Saved",
    isTRUE(input$absq_save_as_new),
    suffix,
    info
  )
  
  removeModal()
})

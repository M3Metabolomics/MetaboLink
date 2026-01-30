  # observeEvent(input$driftMethod, {
  #   if (input$driftMethod == "QC-RFSC (random forrest)") {
  #     hide("dc_qcspan_hide")
  #     hide("dc_degree_hide")
  #     show("dc_ntree_hide")
  #   } else {
  #     hide("dc_ntree_hide")
  #     show("dc_qcspan_hide")
  #     show("dc_degree_hide")
  #   }
  # })
  
observeEvent(input$runDrift, {
  tryCatch({
    # Check if activeFile is NULL
    if (is.null(rv$activeFile)) {
      showNotification("No data - activeFile is NULL", type = "error")
      return()
    }
    
    # Check if sequence exists
    if (is.null(rv$sequence[[rv$activeFile]])) {
      showNotification("No sequence file", type = "error")
      return()
    }
    
    # Get data and sequence
    data <- rv$data[[rv$activeFile]]
    sequence <- rv$sequence[[rv$activeFile]]
    
    # Check if data exists
    if (is.null(data)) {
      showNotification("No data found for active file", type = "error")
      return()
    }
    
    # Check for 'order' column
    if (!"order" %in% colnames(sequence)) {
      showNotification("No 'order' column in sequence, upload sequence with order information", type = "error")
      return()
    }
    
    # Check if order column has any non-NA values
    if (all(is.na(sequence[, 'order']))) {
      showNotification("No order information (all NA), upload sequence with order values", type = "error")
      return()
    }
    
    print(paste("Order column non-NA count:", sum(!is.na(sequence[, 'order']))))
    print(paste("Order column range:", 
                ifelse(all(is.na(sequence[, 'order'])), "All NA", 
                       paste(min(sequence[, 'order'], na.rm = TRUE), 
                             "to", 
                             max(sequence[, 'order'], na.rm = TRUE)))))
    
    # Extract QC data
    qc_indices <- which(sequence[, 1] %in% "QC")
    print(paste("Number of QC samples:", length(qc_indices)))
    
    if (length(qc_indices) == 0) {
      showNotification("No QC samples found in sequence", type = "error")
      return()
    }
    
    dat_qc <- data[, qc_indices, drop = FALSE]
    print(paste("QC data dimensions:", dim(dat_qc)))
    
    # Check if QCs have missing values
    missing_per_qc <- colSums(is.na(dat_qc))
    print(paste("Missing values per QC:", paste(missing_per_qc, collapse = ", ")))
    
    if (any(missing_per_qc > 0)) {
      sendSweetAlert(
        session = session, 
        title = "Error", 
        text = "QCs cannot have missing values. Please impute QC data first.", 
        type = "error"
      )
      return()
    }
    
    # Check if any QC columns are all NA
    if (any(colSums(!is.na(dat_qc)) == 0)) {
      sendSweetAlert(
        session = session, 
        title = "Error", 
        text = "Some QCs have all missing values.", 
        type = "error"
      )
      return()
    }
    
    corrected <- driftCorrection(data, sequence, input$driftMethod, input$driftTrees, input$driftDegree, input$driftQCspan)
    
    rv$tmpData <- corrected
    rv$tmpSequence <- sequence
    updateSelectInput(session, "selectpca1", selected = "Unsaved data", choices = c("Unsaved data", rv$choices))
    output$dttable <- renderDataTable(rv$tmpData, rownames = FALSE, options = list(scrollX = TRUE, scrollY = "700px", pageLength = 20))
    
    print("=== DEBUG DRIFT CORRECTION END ===")
    
  }, error = function(e) {
    showNotification(paste("Error in drift correction:", e$message), type = "error")
    print(paste("Full error:", e))
    print(traceback())  # This will show the call stack
  })
})
  
  observeEvent(input$saveDrift, {
    tryCatch({
      additionalInfo <- paste(
        "Drift correction applied using method: ", input$driftMethod,
        "with ", input$driftTrees, " trees."
      )
      updateDataAndSequence("Drift correct first", input$newFileDrift, "_dc", additionalInfo)
    }, error = function(e) {
      showNotification(paste("Error saving drift correction:", e$message), type = "error")
    })
  })

    observeEvent(input$select_boxplot_1, { #TODO which rv choices -> function
    data <- rv$data[[which(rv$choices %in% input$select_boxplot_1)]]
    sequence <- rv$sequence[[which(rv$choices %in% input$select_boxplot_1)]]
    group <- input$select_boxplot_1_group
    data <- data[sequence[, 1] %in% "Sample" & sequence[, 4] %in% group]
    output$boxplot_1 <- renderPlot({
      boxplot(log2(data), main = input$select_boxplot_1, xlab = "Analyte", ylab = "Intensity")
    })
  }, ignoreInit = TRUE)
  
  observeEvent(input$drift_1, {
    rv$drift_plot_select <- 1
  })
  observeEvent(input$drift_2, {
    rv$drift_plot_select <- 2
  })
  observeEvent(input$drift_3, {
    rv$drift_plot_select <- 3
  })
  
  output$drift_ui <- renderUI({
    box(
      width = NULL,
      if (is.null(rv$activeFile)) {
        p("No data")
      } else if (input$drift_select != "None" && nrow(rv$data[[rv$activeFile]]) != nrow(rv$data[[which(rv$choices %in% input$drift_select)]])) {
        p("Not able to compare the selected datasets")
      } else if (input$drift_select == "None" && rv$drift_plot_select == 2) {
        p("Need dataset to compare with")
      } else if (is.null(input$dt_drift_panel_rows_selected) && rv$drift_plot_select == 1) {
        p("Select feature to plot")
      } else if (rv$drift_plot_select == 1) {
        lapply(seq_along(input$dt_drift_panel_rows_selected), function(i) {
          fluidRow(
            column(6, plotOutput(paste0("driftplotoutput", i), height = 280, width = "100%")),
            column(6, plotOutput(paste0("driftplotoutput2", i), height = 280, width = "100%"))
          )
        })
      } else if (rv$drift_plot_select == 2) {
        fluidRow(column(12, plotOutput("cvscatterplot", height = 600, width = "100%")))
      } else {
        p("Nothing here yet")
      }
    )
  })
  output$boxplot_ui <- renderUI({
    box(
      width = NULL,
      if (is.null(rv$activeFile)) {
        p("No data")
      } else if (is.null(input$dt_boxplot_panel_rows_selected)) {
        p("Select feature to plot")
      } else {
        lapply(seq_along(input$dt_boxplot_panel_rows_selected), function(i) {
          fluidRow(column(12, plotOutput(paste0("boxplotoutput", i), height = 280, width = "100%")))
        })
      }
    )
  })
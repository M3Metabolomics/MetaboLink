observeEvent(input$runImputation, {
  tryCatch({
    
    # Test the condition separately
    cond <- !is.null(rv$activeFile)
    
    # Manually check what need() would return
    if (!cond) {
      print("Condition is FALSE, would trigger 'No data' message")
    } else {
      print("Condition is TRUE, validation passes")
    }
    
    # Try a different approach - use if() instead of validate() temporarily
    if (is.null(rv$activeFile)) {
      showNotification("No data - activeFile is NULL", type = "error")
      return()
    }

    data <- rv$data[[rv$activeFile]]
    sequence <- rv$sequence[[rv$activeFile]]
    
    if (is.null(data) || is.null(sequence)) {
      showNotification("Data or sequence is NULL", type = "error")
      return()
    }
    
    # Check for samples
    if (sum(sequence[, 1] %in% "Sample") < 1) {
      showNotification("Data must have at least one Sample", type = "error")
      return()
    }
    
    # Check for QC - try different columns
    qc_found <- FALSE
    if ("labels" %in% colnames(sequence)) {
      qc_found <- sum(sequence[, "labels"] %in% "QC") >= 1
    } else if (ncol(sequence) >= 1) {
      qc_found <- sum(sequence[, 1] %in% "QC") >= 1
    }
    
    if (!qc_found) {
      showNotification("Data must have at least one QC sample", type = "error")
      return()
    }
    
    print("DEBUG - All checks passed, calling imputation")
    
    imputed <- imputation(data, sequence, input$imputationMethod, input$imputationMinX, input$imp_onlyQC, input$remainingNAs)
    
    rv$tmpData <- imputed
    rv$tmpSequence <- sequence
    updateSelectInput(session, "selectpca1", selected = "Unsaved data", choices = c("Unsaved data", rv$choices))
    output$dttable <- renderDataTable(rv$tmpData, rownames = FALSE, options = list(scrollX = TRUE, scrollY = "700px", pageLength = 20))
    sendSweetAlert(
      title = "Success",
      text = paste0(sum(is.na(rv$data[[rv$activeFile]]) | rv$data[[rv$activeFile]] == 0) - sum(is.na(rv$tmpData) | rv$tmpData == 0), " missing values were imputed."),
      type = "success"
    )
    
    print("=== DEBUG END ===")
    
  }, error = function(e) {
    showNotification(paste("Error in imputation:", e$message), type = "error")
    print(paste("Full error:", e))
    print(traceback())  # This will show the call stack
  })
})
  
  observeEvent(list(input$imputationMethod, input$remainingNAs), {
    tryCatch({
      if (input$imputationMethod == "KNN") {
        hide("imp_minx_hide")
        hide("imp_remaining_hide")
      } else {
        show("imp_remaining_hide")
      }
      if (input$imputationMethod == "Min/X" || input$remainingNAs == "Min/X") {
        show("imp_minx_hide")
      } else {
        hide("imp_minx_hide")
      }
    }, error = function(e) {
      showNotification(paste("Error updating imputation UI:", e$message), type = "error")
    })
  })
  
  observeEvent(input$saveImputation, {
    tryCatch({
      additionalInfo <- paste("Missing values imputation with", input$imputationMethod)
      updateDataAndSequence("Impute first", input$newFileImp, "_imp", additionalInfo)
    }, error = function(e) {
      showNotification(paste("Error saving imputation:", e$message), type = "error")
    })
  })
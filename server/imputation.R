  observeEvent(input$runImputation, {
    tryCatch({
      validate(
        need(!is.null(rv$activeFile), "No data"),
        need(sum(rv$sequence[[rv$activeFile]][, 1] %in% "Sample") >= 1, "Data must have at least one Sample")
      )
      data <- rv$data[[rv$activeFile]]
      sequence <- rv$sequence[[rv$activeFile]]
      validate(need(sum(sequence[, "labels"] %in% "QC") >= 1, "Data must have at least one QC sample"))
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
    }, error = function(e) {
      showNotification(paste("Error in imputation:", e$message), type = "error")
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
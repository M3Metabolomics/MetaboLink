  
  observeEvent(input$blankFiltrate, {
    tryCatch({
      shiny::validate(
        need(!is.null(rv$activeFile), "No data"),
        need("QC" %in% rv$sequence[[rv$activeFile]][, 1], "Data must have at least 1 QC"),
        need("Blank" %in% rv$sequence[[rv$activeFile]][, 1], "Data must have at least 1 Blank"),
        need(sum(rv$sequence[[rv$activeFile]][, 1] %in% "Name") == 1, "Data must have exactly 1 'Name' column")
      )
      sequence <- rv$sequence[[rv$activeFile]]
      data <- rv$data[[rv$activeFile]]
      filtered <- blankFiltration(data, sequence, input$signalStrength, input$keepIS)
      if (input$discardBlank) {
        filtered <- filtered[!sequence[, 1] %in% "Blank"]
        sequence <- sequence[!sequence[, 1] %in% "Blank", ]
      }
      rv$tmpData <- filtered
      rv$tmpSequence <- sequence
      updateSelectInput(session, "selectpca1", selected = "Unsaved data", choices = c("Unsaved data", rv$choices))
      output$dttable <- renderDataTable(rv$tmpData, rownames = FALSE, options = list(scrollX = TRUE, scrollY = "700px", pageLength = 20))
      sendSweetAlert(session, "Success", paste0(nrow(rv$data[[rv$activeFile]]) - nrow(rv$tmpData), " features removed"), type = "success")
    }, error = function(e) {
      showNotification(paste("Error in blank filtration:", e$message), type = "error")
    })
  })
  
  observeEvent(input$saveBF, {
    tryCatch({
      additionalInfo <- paste("Blank filtrated with signal strength above blank =", input$signalStrength)
      updateDataAndSequence("Blank filtrate first", input$newFileBF, paste("_", input$signalStrength, "xb"), additionalInfo)
    }, error = function(e) {
      showNotification(paste("Error saving blank filtration:", e$message), type = "error")
    })
  })
  
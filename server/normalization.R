  #################
  # Normalization #
  #################
  
  observeEvent(input$normalize, {
    if (is.null(rv$activeFile)) {
      showNotification("No data", type = "error")
    } else if(input$normMethod == "QC (PQN)" & sum(rv$sequence[[rv$activeFile]][, 1] %in% "QC") == 0) {
      sendSweetAlert(session = session, title = "Error", text = "No QC samples in dataset.", type = "error")
    } else if(input$normMethod == "Sample amount" & sum(complete.cases(rv$sequence[[rv$activeFile]][, 'amount'])) == 0) {
      sendSweetAlert(session = session, title = "Error", text = "No sample amount information in dataset.", type = "error")
    } else {
      data <- rv$data[[rv$activeFile]]
      sequence <- rv$sequence[[rv$activeFile]]
      qualityControls <- data[, sequence[, 1] %in% "QC"] 
      
      normalizedData <- normalization(data, sequence, qualityControls, input$normMethod)
      
      data[, sequence[, 1] %in% c("QC", "Sample")] <- normalizedData
      
      rv$tmpData <- data
      rv$tmpSequence <- sequence
      
      updateSelectInput(session, "selectpca1", selected = "Unsaved data", 
                        choices = c("Unsaved data", rv$choices))
      output$dttable <- renderDataTable(rv$tmpData, rownames = FALSE, options = 
                                          list(scrollX = TRUE, scrollY = "700px", pageLength = 20))
      sendSweetAlert(title = "Success", text = paste0("Data normalized using ", input$normMethod), type = "success")
    }
  })
  
  observeEvent(input$saveNormalization, {
    additionalInfo <- paste("Normalized with", input$normMethod, " method")
    updateDataAndSequence("Normalize first", input$newFileNorm, "_normalized", additionalInfo)
  })
  
  ##################
  # Transformation #
  ##################
  
  
  observeEvent(input$transform, {
    if (is.null(rv$activeFile)) {
      showNotification("No data", type = "error")
    } else if(input$logTransform == "None" & input$scaling == "None") {
      sendSweetAlert(session = session, title = "Warning", text = "No method selected.", type = "warning")
    } else {
      data <- rv$data[[rv$activeFile]]
      sequence <- rv$sequence[[rv$activeFile]]
      
      cat("Selected transformation method:", input$logTransform, "\n")
      cat("Selected scaling method:", input$scaling, "\n")
      
      transformed <- transformation(data, sequence, input$logTransform, input$scaling)
      
      data[, sequence[, 1] %in% c("QC", "Sample")] <- transformed
      rv$tmpData <- data
      rv$tmpSequence <- sequence
      
      output$dttable <- renderDataTable(rv$tmpData, rownames = FALSE, options = list(scrollX = TRUE, scrollY = "700px", pageLength = 20))
      sendSweetAlert(session, title = "Success", text = "Data transformed.", type = "success")
    }
  })
  
  observeEvent(input$saveTransform, {
    additionalInfo <- paste(
      "Transformed with log transformation: ", input$logTransform,
      " and scaling: ", input$scaling
    )
    updateDataAndSequence("Transform first", input$newFileTransform, "_transformed", additionalInfo)
  })
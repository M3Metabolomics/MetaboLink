observeEvent(input$normalizeIS, {
  tryCatch({
    
    # Check if activeFile is NULL
    if (is.null(rv$activeFile)) {
      showNotification("No data - activeFile is NULL", type = "error")
      return()
    }
    
    # Get data and sequence
    data <- rv$data[[rv$activeFile]]
    sequence <- rv$sequence[[rv$activeFile]]
    
    # Check if data and sequence exist
    if (is.null(data) || is.null(sequence)) {
      showNotification("Data or sequence is NULL", type = "error")
      return()
    }
    
    # Check for exactly 1 'Name' column in the first column of sequence
    name_count <- sum(sequence[, 1] %in% "Name")
    print(paste("'Name' count in first column:", name_count))
    
    if (name_count != 1) {
      showNotification("Data must have exactly 1 'Name' column", type = "error")
      return()
    }
    
    # Check if internal standards are selected
    if (is.null(input$isChoose)) {
      showNotification("No internal standards selected", type = "error")
      return()
    }
    
    # Call normalization function
    normalized <- normalizationIS(data, sequence, input$isChoose, input$isMethod, input$normalizeQC)
    
    if (is.null(normalized)) {
      sendSweetAlert(session, "Error", "Internal standard normalization failed due to missing values in IS.", type = "error")
    } else {
      # Add isnorm column to sequence
      isColumn <- c("-", rep(NA, ncol(sequence) - 1))
      sequence <- rbind(sequence, isColumn)
      rownames(sequence)[nrow(sequence)] <- "isnorm"
      
      rv$tmpData <- normalized
      rv$tmpSequence <- sequence
      
      sendSweetAlert(
        session, 
        title = "Success", 
        text = paste0("Internal standards normalized with ", input$isMethod, " method"), 
        type = "success"
      )
      
      updateSelectInput(session, "selectpca1", selected = "Unsaved data", choices = c("Unsaved data", rv$choices))
      output$dttable <- renderDataTable(rv$tmpData, rownames = FALSE, options = list(scrollX = TRUE, scrollY = "700px", pageLength = 20))
    }
    
    print("=== DEBUG IS NORMALIZATION END ===")
    
  }, error = function(e) {
    showNotification(paste("Error in normalization:", e$message), type = "error")
    print(paste("Full error:", e))
    print(traceback())  # This will show the call stack
  })
})
  
  observeEvent(input$saveIS, {
    tryCatch({
      additionalInfo <- paste("Internal standards normalized with", input$isMethod, "method")
      updateDataAndSequence("IS normalize first", input$newFileIS, "_is", additionalInfo)
    }, error = function(e) {
      showNotification(paste("Error saving IS normalization:", e$message), type = "error")
    })
  })
  
  observeEvent(input$removeIS, {
    tryCatch({
      validate(need(!is.null(rv$activeFile), "No data"))
      data <- rv$data[[rv$activeFile]]
      sequence <- rv$sequence[[rv$activeFile]]
      toRemove <- data[sequence[, 1] %in% "Name"]
      data <- data[!grepl("\\(IS\\)", toupper(toRemove[ , 1])), ]
      rv$data[[rv$activeFile]] <- data
      updateCheckboxGroupInput(session, "isChoose", choices = character(0), selected = NULL)
    }, error = function(e) {
      showNotification(paste("Error removing IS:", e$message), type = "error")
    })
  })
  
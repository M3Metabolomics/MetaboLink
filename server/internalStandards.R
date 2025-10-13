  observeEvent(input$normalizeIS, {
    tryCatch({
      validate(
        need(!is.null(rv$activeFile), "No data"),
        need(sum(rv$sequence[[rv$activeFile]][, 1] %in% "Name") == 1, "Data must have exactly 1 'Name' column"),
        need(!is.null(input$isChoose), "No internal standards selected")
      )
      sequence <- rv$sequence[[rv$activeFile]]
      data <- rv$data[[rv$activeFile]]
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
        sendSweetAlert(session, title = "Success", text = paste0("Internal standards normalized with ", input$isMethod, " method"), type = "success")
        updateSelectInput(session, "selectpca1", selected = "Unsaved data", choices = c("Unsaved data", rv$choices))
        output$dttable <- renderDataTable(rv$tmpData, rownames = FALSE, options = list(scrollX = TRUE, scrollY = "700px", pageLength = 20))
      }
    }, error = function(e) {
      showNotification(paste("Error in normalization:", e$message), type = "error")
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
  
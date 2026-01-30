  
observeEvent(input$runFilterNA, {
  # SIMPLIFIED VERSION - Remove all complex checks
  
  tryCatch({
    # 1. Basic validation
    if (is.null(rv$activeFile)) {
      sendSweetAlert(session, "Error!", "No data loaded.", type = "error")
      return()
    }
    
    # 2. Get data and sequence
    current_data <- rv$data[[rv$activeFile]]
    sequence <- rv$sequence[[rv$activeFile]]
    method <- input$filterNAmethod
    
    # 3. Simple method check
    if (is.null(method) || length(method) == 0) {
      sendSweetAlert(session, "Error!", "Please select at least one method.", type = "error")
      return()
    }
    
    # 4. Check for group method (simplified)
    if ("in group" %in% method) {
      if (is.null(sequence) || ncol(sequence) < 4 || all(is.na(sequence[, 4]))) {
        sendSweetAlert(session, "Error!", "Valid group information is required.", type = "error")
        return()
      }
    }
    
    # 5. Run filtration
    mvf_dat <- cutoffrm(current_data, sequence, input$cutoffNAs, method)
    
    # 6. Store results
    rv$tmpData <- mvf_dat
    rv$tmpSequence <- sequence
    
    # 7. Update UI
    updateSelectInput(session, "selectpca1", selected = "Unsaved data", 
                      choices = c("Unsaved data", rv$choices))
    
    output$dttable <- renderDataTable(
      rv$tmpData, 
      rownames = FALSE, 
      options = list(scrollX = TRUE, scrollY = "700px", pageLength = 20)
    )
    
    sendSweetAlert(
      title = "Success",
      text = paste0(nrow(current_data) - nrow(rv$tmpData), " feature(s) removed"),
      type = "success"
    )
    
  }, error = function(e) {
    # Detailed error
    err_msg <- paste("Error:", e$message)
    cat("FINAL ERROR:", err_msg, "\n")
    showNotification(err_msg, type = "error")
  })
})
  
  observeEvent(input$saveFilterNA, {
    tryCatch({
      additionalInfo <- paste(
        "Missing value filtration using",
        input$cutoffNAs,
        "% as threshold and method -",
        paste(input$filterNAmethod, collapse=", ")
      )
      updateDataAndSequence("Filtrate first", input$mvf_newsave, "_mvr", additionalInfo)
    }, error = function(e) {
      showNotification(paste("Error saving missing value filtration:", e$message), type = "error")
    })
  })
  
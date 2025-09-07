  
  observeEvent(input$runFilterNA, {
    tryCatch({
      validate(need(!is.null(rv$activeFile), "No data"))
      sequence <- rv$sequence[[rv$activeFile]]
      method <- input$filterNAmethod
      if (("in group" %in% method) & !any(complete.cases(sequence[, 4]))) {
        sendSweetAlert(session, "Error!", "Group information needed.", type = "error")
      } else if (is.null(method)) {
        sendSweetAlert(session, "Error!", "No method selected.", type = "error")
      } else {
        mvf_dat <- cutoffrm(rv$data[[rv$activeFile]], sequence, input$cutoffNAs, method)
        rv$tmpData <- mvf_dat
        rv$tmpSequence <- sequence
        updateSelectInput(session, "selectpca1", selected = "Unsaved data", choices = c("Unsaved data", rv$choices))
        output$dttable <- renderDataTable(rv$tmpData, rownames = FALSE, options = list(scrollX = TRUE, scrollY = "700px", pageLength = 20))
        sendSweetAlert(
          title = "Success",
          text = paste0(nrow(rv$data[[rv$activeFile]]) - nrow(rv$tmpData), " feature(s) removed"),
          type = "success"
        )
      }
    }, error = function(e) {
      showNotification(paste("Error in missing value filtration:", e$message), type = "error")
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
  
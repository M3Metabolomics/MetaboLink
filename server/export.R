observeEvent(input$export_xml_list, {
  tryCatch({
    # Check if there is anything selected
    if (length(input$export_xml_list) == 0) {
      showNotification("No data selected for export", type = "error")
      return()
    }
    # Create the download handler
    output$export_xml <- downloadHandler(
      filename = function() {
        # Just use the first selected item for the filename
        paste0(input$export_xml_list[1], ".xlsx")
      },
      content = function(file) {
        # Get the indices of selected data
        indices <- which(rv$choices %in% input$export_xml_list)
        
        # Subset the data
        export_data <- rv$data[indices]
        
        # Write to Excel
        write_xlsx(export_data, file)
      }
    )
    
  }, error = function(e) {
    showNotification(paste("Error exporting:", e$message), type = "error")
  })
})
#########################
# Volcano Results Export #
#########################

# Render volcano export buttons in a separate UI output
output$volcano_export_buttons <- renderUI({
  if (!is.null(rv$volcano_export_names) && length(rv$volcano_export_names) > 0) {
    tagList(
      h5("Volcano Results:"),
      lapply(seq_along(rv$volcano_export_names), function(i) {
        downloadButton(paste0("dwn_volcano_", i), 
                       paste0(rv$volcano_export_names[i], ".xlsx"),
                       class = "btn-sm btn-info",
                       style = "margin-bottom: 5px;")
      })
    )
  } else {
    return(NULL)
  }
})

# Create download handlers for volcano results
observe({
  if (!is.null(rv$volcano_export_names) && length(rv$volcano_export_names) > 0) {
    for (i in seq_along(rv$volcano_export_names)) {
      local({
        idx <- i
        name <- rv$volcano_export_names[idx]
        output[[paste0("dwn_volcano_", idx)]] <- downloadHandler(
          filename = function() {
            paste0(name, ".xlsx")
          },
          content = function(file) {
            write_xlsx(rv$volcano_exports[[name]], file)
          }
        )
      })
    }
  }
})
#########################
# PolySTest and VSClust #
#########################

observeEvent(input$export_polystest, {
  tryCatch({
    validate(
      need(!is.null(rv$activeFile), "No data loaded")
    )
    sequence <- rv$sequence[[rv$activeFile]]
    tdata <- rv$data[[rv$activeFile]][, sequence[, 1] %in% c("Name",  "Sample")]
    groups <- c(input$group1_polystest, input$group2_polystest)
    time <- c(input$time1_polystest, input$time2_polystest)
    selected <- selectPolySTest(tdata, sequence, groups, time)
    PolySTestMessage <- prepareMessage2(selected$selected, selected$selected_sequence, time)
    js$send_message(url="http://computproteomics.bmb.sdu.dk:443/app_direct/PolySTest/", 
                    dat=PolySTestMessage, tool="PolySTest")
  }, error = function(e) {
    showNotification(paste("Error exporting PolySTest:", e$message), type = "error")
  })
})

observeEvent(input$send_polystest, {
  tryCatch({
    validate(
      need(!is.null(rv$activeFile), "No data loaded")
    )
    sequence <- rv$sequence[[rv$activeFile]]
    tdata <- rv$data[[rv$activeFile]][, sequence[, 1] %in% c("Name",  "Sample")]
    tseq <- sequence[sequence[, 1] %in% c("Name",  "Sample"), ]
    time <- complete.cases(tseq[, 5])
    if(any(complete.cases(tseq[, 5]))) {
      time <- unique(tseq[complete.cases(tseq[, 5]), 5])
    } else {
      time <- c("")
    }
    PolySTestMessage <- prepareMessage2(tdata, tseq, time)
    js$send_message(url="http://computproteomics.bmb.sdu.dk:443/app_direct/PolySTest/", 
                    dat=PolySTestMessage, tool="PolySTest")
  }, error = function(e) {
    showNotification(paste("Error sending PolySTest:", e$message), type = "error")
  })
})

observeEvent(input$send_vsclust, {
  tryCatch({
    validate(
      need(!is.null(rv$activeFile), "No data loaded")
    )
    sequence <- rv$sequence[[rv$activeFile]]
    tdata <- rv$data[[rv$activeFile]][, sequence[, 1] %in% c("Name",  "Sample")]
    tseq <- sequence[sequence[, 1] %in% c("Name",  "Sample"), ]
    VSClustMessage <- prepareMessage2(tdata, tseq)
    js$send_message(url="http://computproteomics.bmb.sdu.dk/app_direct/VSClust/",
                    dat=VSClustMessage, tool="VSClust")
  }, error = function(e) {
    showNotification(paste("Error sending VSClust:", e$message), type = "error")
  })
})


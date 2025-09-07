  observeEvent(input$export_xml_list, {
    tryCatch({
      validate(
        need(length(input$export_xml_list) > 0, "No data selected for export")
      )
      output$export_xml <- downloadHandler(
        filename = function() {
          paste0(names(rv$data[rv$choices %in% input$export_xml_list])[1], ".xlsx")
        },
        content = function(file) {
          write_xlsx(rv$data[rv$choices %in% input$export_xml_list], file)
        }
      )
    }, error = function(e) {
      showNotification(paste("Error exporting XML:", e$message), type = "error")
    })
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
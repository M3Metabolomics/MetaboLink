  
  observeEvent(input$inputFile, { # Ensure file is uploaded before pressing Upload button
    inputFile <<- read.csv(input$inputFile$datapath,
                           header = 1,
                           stringsAsFactors = F,
                           check.names = FALSE)
    enable("upload")
  })
  
  observeEvent(input$upload, {
    shinyCatch({
      if(input$fileType == "Samples in rows") {
        #TODO
        inputFile <- t(inputFile)
        inputFile <- data.frame(inputFile, stringsAsFactors = FALSE)
      }
      
      # Look for column "Name" or "name" and make it the first column
      if("Name" %in% colnames(inputFile)) {
        inputFile <- inputFile[, c("Name", setdiff(colnames(inputFile), "Name"))]
      } else if("name" %in% colnames(inputFile)) {
        inputFile <- inputFile[, c("name", setdiff(colnames(inputFile), "name"))]
      }
      
      # Remove the special charactor special ± character (Unicode U+00B1) from data[,1]
      inputFile[,1] <- iconv(inputFile[,1], "WINDOWS-1252", "UTF-8", sub = "")
      inputFile[,1] <- gsub("\\(±\\)", "", inputFile[,1])
    
    },
    blocking_level = 'message'
    )
    if(any(duplicated(names(inputFile)))) {
      sendSweetAlert(session,
                     title = "Error",
                     text = paste("Duplicate columns found."),
                     type = "error")
    } else {
      labels <- identifyLabels(inputFile)
      initializeVariables()
      rv$sequence[[length(rv$sequence) + 1]] <- data.frame(labels,
                                                           batch = NA,
                                                           order = NA,
                                                           group = NA,
                                                           time = NA,
                                                           paired = NA,
                                                           amount = NA)
      rv$data[[length(rv$data) + 1]] <- inputFile
      names(rv$data)[length(rv$data)] <- substr(input$inputFile$name, 1, nchar(input$inputFile$name) - 4)
      
      rv$choices <- paste(seq_along(rv$data), ": ", names(rv$data))
      rv$activeFile <- length(rv$data)
      updateTabItems(session, "tabs", selected = "Datainput")
      show("buttons")
    }
  })
  
  observeEvent(input$inputSequence, {
    shinyCatch({
      inputSequence <- read.csv(input$inputSequence$datapath, header = 1, stringsAsFactors = FALSE)
      colnames(inputSequence) <- tolower(colnames(inputSequence))
      inputSequence <- checkSequence(inputSequence)
    },
    blocking_level = 'message'
    )
    sequence <- rv$sequence[[rv$activeFile]]
    labeledSequence <- data.frame("sample" = row.names(sequence), sequence)
    inputSequence["sample"] <- lapply(inputSequence["sample"], as.character)
    sequence <- left_join(labeledSequence[, 1:2], inputSequence, by = "sample")
    row.names(sequence) <- sequence[, 1]
    sequence <- sequence[, -1]
    rv$sequence[[rv$activeFile]] <- sequence
    
    if(any(grepl("[^[:alnum:]_]", sequence$group))) {
      showModal(
        modalDialog(
          title = "Invalid group names", size = "m", easyClose = TRUE,
          footer = list(actionButton("group_name_format", "Format names"), modalButton("Dismiss")),
          fluidRow(
            column(12, p("Invalid group names found. Group names must be alphanumeric and not include spaces. The use of ´_´ has been allowed."))
          )
        )
      )
    }
  })
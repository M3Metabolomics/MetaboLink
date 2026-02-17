  #TODO
  observeEvent(input$run_pca1, {
    if (!is.null(rv$activeFile)) { 
      if (input$selectpca1 == "Unsaved data") {
        data <- rv$tmpData       # Set data to the temporary data
        seq <- rv$tmpSequence    # Set sequence to the temporary sequence
      } else { 
        selectchoices <- paste(seq_along(rv$data), ": ", names(rv$data)) # Get the selected dataset
        sd <- which(rv$choices %in% input$selectpca1) # Get the index of the selected dataset
        data <- rv$data[[sd]]    # Set data to the selected dataset
        seq <- rv$sequence[[sd]] # Set sequence to the selected sequence
      }
      
      if ("Sample"  %in% seq[, "labels"]) { # Check if the sequence file contains a "Sample" column
        if (any(seq[, "labels"] %in% "QC")) { # Check if the sequence file contains a "QC" column
          seq[seq[, "labels"] %in% "QC", "group"] <- "QC" # Set the "QC" column to "QC"
        } else {
          cat("No 'QC' labels found in the sequence.\n")
        }
        
        data_subset <- data[seq[, "labels"] %in% c("Sample", "QC")] # Get the data for the samples and QC
        
        
        # Check if Name column exists
        if (!"Name" %in% colnames(data)) {
          showModal(
            modalDialog(
              title = "Missing Name Column", 
              size = "m", 
              easyClose = TRUE,
              footer = list(
                actionButton("assign_names_pca1", "Create Name column", 
                             icon = icon("plus"), 
                             style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                modalButton("Dismiss", icon = icon("times"))
              ),
              fluidRow(
                column(12, 
                       p("The dataset does not have a 'Name' column."),
                       p("Feature names are required for PCA visualization and interpretation."),
                       p("Would you like to create a Name column with placeholder names?")
                )
              )
            )
          )
          return()
        }
        
        # Check if Name column has empty or NA values
        if (any(is.na(data[, "Name"]) | data[, "Name"] == "")) {
          showModal(
            modalDialog(
              title = "Missing Feature Names", 
              size = "m", 
              easyClose = TRUE,
              footer = list(
                actionButton("assign_names_pca1", "Assign names", 
                             icon = icon("pencil"), 
                             style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                modalButton("Dismiss", icon = icon("times"))
              ),
              fluidRow(
                column(12, 
                       p("The dataset has missing or empty names in the 'Name' column."),
                       p("Feature names are required for PCA visualization and interpretation.")
                )
              ),
              br(),
              fluidRow(
                column(6, 
                       wellPanel(
                         h5("Assign Names", style = "color: #337ab7;"),
                         p("Generate placeholder names", style = "font-size: 12px;"),
                         p("(e.g., Feature1, Feature2, ...)", style = "font-size: 12px; font-style: italic;")
                       )
                ),
                column(6, 
                       wellPanel(
                         h5("Dismiss", style = "color: #777;"),
                         p("Cancel PCA generation", style = "font-size: 12px;"),
                         p("and fix names manually", style = "font-size: 12px; font-style: italic;")
                       )
                )
              )
            )
          )
          return()
        }
        rownames(data_subset) <- make.unique(as.character(data[, "Name"])) # Make the rownames unique
        
        seq_subset <- seq[seq[, "labels"] %in% c("Sample", "QC"), ] # Get the sequence for the samples and QC
        
        # Perform PCA once and save the results to pca_result
        pca_result <- pcaplot(data_subset, seq_subset, input$pca1_islog)
        
        message("PCA results saved.")
        # Generate a unique name for the PCA result based on the dataset name
        if (input$selectpca1 == "Unsaved data") {
          dataset_name <- "UnsavedData"  # or any other name you prefer for unsaved data
        } else {
          dataset_name <- names(rv$data)[sd]
        }
        pca_name <- paste0(dataset_name, "_pca")
        pc_name <- paste0(dataset_name, "_PC")
        
        # Check if the PCA name already exists in rv$pca_results
        if (!(pca_name %in% names(rv$pca_results))) {
          # If the name does not exist, save the PCA and PC results
          rv$pca_results[[pca_name]] <- list(pca_df = pca_result$pca_df,
                                             PC_df = pca_result$PC_df)
        }
        
        output$plotpca1 <- renderPlotly({
          pca_result$pca_plotly
        })
        
        output$plotscree1 <- renderPlotly({
          pca_result$scree_plotly
        })
        
        message(sample(quotes, 1))
        
        if (sum(seq[, 1] %in% "QC") > 0) {
          qccv <- paste0("CV in QC samples: ", round(cvmean(data[seq[, 1] %in% "QC"]), 2), "</br>")
        } else {
          qccv <- "No QC in dataset </br>"
        }
        sclass <- seq[seq[, "labels"] %in% c("Sample", "QC"), ][, "group"] # Get the class of the samples and QC
        sclass <- sclass[sclass != "QC"]
        if (sum(!is.na(sclass)) > 0) {
          classcv <- sapply(sort(unique(sclass)), function(x) {
            round(cvmean(data_subset[, sclass %in% x]), 2)
          })
          classcv <- sapply(seq_along(classcv), function(x) {
            paste0("CV in group ", sort(unique(sclass))[x], ": ", classcv[x], "</br>")
          })
        } else {
          classcv <- NULL
        }
        text <- c(qccv, classcv)
        output$pca1Details <- renderUI({
          HTML(text)
        })
      }
    }
  })
  
  # Observer for assigning names in PCA1
  observeEvent(input$assign_names_pca1, {
    # Remove the modal
    removeModal()
    
    tryCatch({
      # Get current data based on selection
      if (input$selectpca1 == "Unsaved data") {
        data <- rv$tmpData
        seq <- rv$tmpSequence
      } else {
        sd <- which(rv$choices %in% input$selectpca1)
        data <- rv$data[[sd]]
        seq <- rv$sequence[[sd]]
      }
      
      # Check if Name column exists
      if (!"Name" %in% colnames(data)) {
        # Create Name column with placeholder names
        data$Name <- paste0("Feature", 1:nrow(data))
        # Reorder columns to put Name first
        data <- data[, c("Name", setdiff(colnames(data), "Name"))]
      } else {
        # Assign placeholder names where Name is empty or NA
        missing_name_idx <- which(is.na(data[, "Name"]) | data[, "Name"] == "")
        
        # Create placeholder names (Feature + row number)
        for (i in seq_along(missing_name_idx)) {
          data[missing_name_idx[i], "Name"] <- paste0("Feature", missing_name_idx[i])
        }
      }
      
      # Update the data
      if (input$selectpca1 == "Unsaved data") {
        rv$tmpData <- data
        rv$tmpSequence <- seq
      } else {
        sd <- which(rv$choices %in% input$selectpca1)
        rv$data[[sd]] <- data
        rv$sequence[[sd]] <- seq
      }
      
      # Show success message
      sendSweetAlert(
        session, 
        "Success", 
        ifelse(! "Name" %in% colnames(data),
               "Name column created with placeholder names.",
               paste0("Assigned placeholder names to ", length(missing_name_idx), " features.")),
        type = "success"
      )
      
      # Automatically re-run PCA after assigning names
      # This will trigger the PCA calculation again
      shinyjs::click("run_pca1")
      
    }, error = function(e) {
      showNotification(paste("Error assigning names:", e$message), type = "error")
    })
  })
  
  observeEvent(input$run_pca2, {
    selectchoices <- paste(seq_along(rv$data), ": ", names(rv$data))
    sd <- which(rv$choices %in% input$selectpca2)
    if ("Sample" %in% rv$sequence[[sd]][, 1]) {
      data <- rv$data[[sd]]
      seq <- rv$sequence[[sd]]
      shinyCatch(
        seq[seq[, 1] %in% "QC", ][, 4] <- "QC",
        blocking_level = 'message',
        shiny = FALSE
      )
      
      data_subset <- data[seq[, "labels"] %in% c("Sample", "QC")] # Get the data for the samples and QC
      # Check if Name column exists
      if (!"Name" %in% colnames(data)) {
        showModal(
          modalDialog(
            title = "Missing Name Column", 
            size = "m", 
            easyClose = TRUE,
            footer = list(
              actionButton("assign_names_pca2", "Create Name column", 
                           icon = icon("plus"), 
                           style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
              modalButton("Dismiss", icon = icon("times"))
            ),
            fluidRow(
              column(12, 
                     p("The dataset does not have a 'Name' column."),
                     p("Feature names are required for PCA visualization and interpretation."),
                     p("Would you like to create a Name column with placeholder names?")
              )
            )
          )
        )
        return()
      }
      
      # Check if Name column has empty or NA values
      if (any(is.na(data[, "Name"]) | data[, "Name"] == "")) {
        showModal(
          modalDialog(
            title = "Missing Feature Names", 
            size = "m", 
            easyClose = TRUE,
            footer = list(
              actionButton("assign_names_pca2", "Assign names", 
                           icon = icon("pencil"), 
                           style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
              modalButton("Dismiss", icon = icon("times"))
            ),
            fluidRow(
              column(12, 
                     p("The dataset has missing or empty names in the 'Name' column."),
                     p("Feature names are required for PCA visualization and interpretation.")
              )
            ),
            br(),
            fluidRow(
              column(6, 
                     wellPanel(
                       h5("Assign Names", style = "color: #337ab7;"),
                       p("Generate placeholder names", style = "font-size: 12px;"),
                       p("(e.g., Feature1, Feature2, ...)", style = "font-size: 12px; font-style: italic;")
                     )
              ),
              column(6, 
                     wellPanel(
                       h5("Dismiss", style = "color: #777;"),
                       p("Cancel PCA generation", style = "font-size: 12px;"),
                       p("and fix names manually", style = "font-size: 12px; font-style: italic;")
                     )
              )
            )
          )
        )
        return()
      }
      rownames(data_subset) <- make.unique(as.character(data[, "Name"]))
      
      seq_subset <- seq[seq[, "labels"] %in% c("Sample", "QC"), ] # Get the sequence for the samples and QC
      
      # Save the PCA results to pca_results
      pca_result <- pcaplot(data_subset, seq_subset, input$pca2_islog)
      
      message("PCA results saved.")
      # Generate a unique name for the PCA result based on the dataset name
      if (input$selectpca1 == "Unsaved data") {
        dataset_name <- "UnsavedData"  # or any other name you prefer for unsaved data
      } else {
        dataset_name <- names(rv$data)[sd]
      }
      pca_name <- paste0(dataset_name, "_pca")
      pc_name <- paste0(dataset_name, "_PC")
      
      # Check if the PCA name already exists in rv$pca_results
      if (!(pca_name %in% names(rv$pca_results))) {
        # If the name does not exist, save the PCA and PC results
        pca_result <- pcaplot(data_subset, seq_subset, input$pca2_islog)  # Perform PCA
        
        # Save the PCA and PC results as a named list for each PCA result
        rv$pca_results[[pca_name]] <- list(pca_df = pca_result$pca_df,
                                           PC_df = pca_result$PC_df)
      }
      
      # Debugging to show that rv$results is updated
      # cat("PCA results saved as:", pca_name, "\n")
      # cat("PCA results dimensions:", dim(rv$pca_results[[pca_name]]), "\n")
      # print(str(rv$pca_results[[pca_name]]))
      
      output$plotpca2 <- renderPlotly({
        pca_result$pca_plotly
      })
      
      output$plotscree2 <- renderPlotly({
        pca_result$scree_plotly
      })
      
      if (sum(seq$labels %in% "QC") > 0) {
        qccv <- paste0("CV in QC samples: ", round(cvmean(data[seq[, 1] %in% "QC"]), 2), "</br>")
      } else {
        qccv <- "No QC in dataset </br>"
      }
      sclass <- seq[seq[, 1] %in% c("Sample", "QC"), ][, 4]
      sclass <- sclass[sclass != "QC"]
      if (sum(!is.na(sclass)) > 0) {
        classcv <- sapply(sort(unique(sclass)), function(x) {
          round(cvmean(data_subset[sclass %in% x]), 2)
        })
        classcv <- sapply(seq_along(classcv), function(x) {
          paste0("CV in group ", sort(unique(sclass))[x], ": ", classcv[x], "</br>")
        })
      } else {
        classcv <- NULL
      }
      text <- c(qccv, classcv)
      output$pca2Details <- renderUI({
        HTML(text)
      })
    }
  })
  
  # Observer for assigning names in PCA2
  observeEvent(input$assign_names_pca2, {
    # Remove the modal
    removeModal()
    
    tryCatch({
      # Get current data based on selection
      sd <- which(rv$choices %in% input$selectpca2)
      data <- rv$data[[sd]]
      seq <- rv$sequence[[sd]]
      
      # Check if Name column exists
      if (!"Name" %in% colnames(data)) {
        # Create Name column with placeholder names
        data$Name <- paste0("Feature", 1:nrow(data))
        # Reorder columns to put Name first
        data <- data[, c("Name", setdiff(colnames(data), "Name"))]
      } else {
        # Assign placeholder names where Name is empty or NA
        missing_name_idx <- which(is.na(data[, "Name"]) | data[, "Name"] == "")
        
        # Create placeholder names (Feature + row number)
        for (i in seq_along(missing_name_idx)) {
          data[missing_name_idx[i], "Name"] <- paste0("Feature", missing_name_idx[i])
        }
      }
      
      # Update the data
      rv$data[[sd]] <- data
      rv$sequence[[sd]] <- seq
      
      # Show success message
      sendSweetAlert(
        session, 
        "Success", 
        ifelse(! "Name" %in% colnames(data),
               "Name column created with placeholder names.",
               paste0("Assigned placeholder names to ", length(missing_name_idx), " features.")),
        type = "success"
      )
      
      # Automatically re-run PCA after assigning names
      shinyjs::click("run_pca2")
      
    }, error = function(e) {
      showNotification(paste("Error assigning names:", e$message), type = "error")
    })
  })
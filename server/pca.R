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
        
        if (any(is.na(data[, "Name"]) | data[, "Name"] == "")) {
          sendSweetAlert(session, "Error",
                         "No names in Name column.
                         Make sure features has names ;)",
                         type = "error")
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
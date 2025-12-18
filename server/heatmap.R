  # Generate Heatmap
  output$group_selection_ui_heatmap <- renderUI({
    
    if (!is.null(rv$activeFile)) {
      if (input$select_heatmap_data == "Unsaved data") {
        data <- rv$tmpData  # Use the temporary data
        seq <- rv$tmpSequence  # Use the temporary sequence
      } else {
        # Get the index of the selected dataset
        sd <- which(rv$choices %in% input$select_heatmap_data)
        data <- rv$data[[sd]]  # Retrieve the selected dataset
        seq <- rv$sequence[[sd]]  # Retrieve the selected sequence
      }
      
      if (input$select_groups_heatmap) {  # Only render if the checkbox is checked
        selectInput(
          "selected_groups_heatmap", 
          "Select Groups:", 
          choices = seq$group, 
          selected = seq$group[1],  # Default to the first group
          multiple = TRUE,          # Allow multiple selections
          width = "100%"
        )
      }
    }
  })
  output$grouping_column_ui <- renderUI({
    if (!is.null(rv$activeFile)) {
      if (input$select_heatmap_data == "Unsaved data") {
        data <- rv$tmpData  # Use the temporary data
        seq <- rv$tmpSequence  # Use the temporary sequence
      } else {
        # Get the index of the selected dataset
        sd <- which(rv$choices %in% input$select_heatmap_data)
        data <- rv$data[[sd]]  # Retrieve the selected dataset
        seq <- rv$sequence[[sd]]  # Retrieve the selected sequence
      }
      
      # seq <- seq[!seq[, "labels"] %in% c("Sample", "QC"), ]  # Restrict rows to "Sample" and "QC"
      # data <- data[, rownames(seq), drop = FALSE]  # Use row names of seq to filter columns
      
      columns <- colnames(data)
      
      if ("super_class" %in% columns) {
        default_column <- "super_class"
      } else {
        default_column <- columns[1]
      }
      
      if (input$enable_grouping_heatmap) {  # Only render if the checkbox is checked
        selectInput(
          inputId = "group_column_heatmap",
          label = "Select grouping column",
          choices = columns,
          selected = default_column, 
          width = "100%"
        )
      }
    }
  })
  
  # Define the reactive value at the top of the server so it persists
  savedDatasetNameHeatmap <- reactiveVal("My Heatmap")
  # Observe the input for the heatmap title and update the reactive value
  observe({
    savedDatasetNameHeatmap(input$heatmap_title)
    
    output$displayName <- renderText({
      paste("Current Heatmap Title:", savedDatasetNameHeatmap())
    })
    
  })
  
  observeEvent({input$select_heatmap_data}, {
    if (!is.null(rv$activeFile)) {
      if (input$select_heatmap_data == "Unsaved data") {
        data <- rv$tmpData  # Use the temporary data
        seq <- rv$tmpSequence  # Use the temporary sequence
      } else {
        # Get the index of the selected dataset
        sd <- which(rv$choices %in% input$select_heatmap_data)
        data <- rv$data[[sd]]  # Retrieve the selected dataset
        seq <- rv$sequence[[sd]]  # Retrieve the selected sequence
      }
      
      # TODO crazes the example file
      # seq <- seq[!seq[, "labels"] %in% c("Sample", "QC"), ] 
      # data_sub <- data[, rownames(seq), drop = FALSE]  # Use row names of seq to filter columns
      
      # # Extract column names from the selected dataset
      data_colnames <- colnames(data) # Substitude with data_sub
      columns <- c("heatmap_labels")
      for (column in columns) {
        # Update the 'identifier_column' select input with the new choices
        updateSelectInput(session, column, choices = data_colnames)
      }
    }
  })
  
  observeEvent(input$run_heatmap, {
    # Ensure a dataset is selected
    req(input$select_heatmap_data,
        input$heatmap_labels)
    
    if (!is.null(rv$activeFile)) {
      if (input$select_heatmap_data == "Unsaved data") {
        data <- rv$tmpData  # Use the temporary data
        seq <- rv$tmpSequence  # Use the temporary sequence
      } else {
        # Get the index of the selected dataset
        sd <- which(rv$choices %in% input$select_heatmap_data)
        data <- rv$data[[sd]]  # Retrieve the selected dataset
        seq <- rv$sequence[[sd]]  # Retrieve the selected sequence
        # dataset_name <- names(rv$data)[sd]  # Retrieve dataset name
      }
      
      # Subset data for "Sample" labels
      seq_subset <- seq[seq[, "labels"] %in% c("Sample", 2), ]  # Restrict to "Sample" rows
      data_subset <- data[, rownames(seq_subset), drop = FALSE]  # Use row names of seq_subset to filter columns
      
      # Check group selection
      if (input$select_groups_heatmap) {
        if (is.null(input$selected_groups_heatmap) || length(input$selected_groups_heatmap) < 2) {
          showNotification("Please select at least two groups for the heatmap.", type = "error")
          return()  # Stop execution
        }
        # Filter seq_subset and data_subset by selected groups
        selected_groups_heatmap <- input$selected_groups_heatmap
        seq_subset <- seq_subset[seq_subset$group %in% selected_groups_heatmap, ]
        data_subset <- data[, rownames(seq_subset), drop = FALSE]  # Subset columns by rownames of seq_subset
      }
      
      
      enable_groups <- input$enable_grouping_heatmap
      groups <- input$group_column_heatmap
      show_column_names <- input$show_column_names
      show_row_names <- input$show_row_names
      cluster_rows <- input$cluster_rows
      show_row_dend <- input$show_row_dend
      labels <- input$heatmap_labels
      clustering_distance_rows <- input$clustering_distance_rows
      clustering_method_rows <- input$clustering_method_rows
      islog <- input$heatmap_islog
      
      message(paste0("Heatmap labels column: ", labels))
      
      message(paste0("Enable grouping: ", enable_groups))
      if (enable_groups) {
        message(paste0("Grouping column selected: ", groups))
      }
      
      
      selected_labels <- as.character(data[[labels]])
      fallback <- if ("Name" %in% colnames(data)) {
        as.character(data[["Name"]])
      } else if ("name" %in% colnames(data)) {
        as.character(data[["name"]])
      } else {
        NULL
      }
      if (is.null(fallback)) {
        showNotification("No fallback column ('Name' or 'name') available.", type = "error")
        return()
      }
      missing <- is.na(selected_labels) | selected_labels == ""
      selected_labels[missing] <- fallback[missing]
      rownames(data_subset) <- make.unique(selected_labels)
      rownames(data) <- make.unique(selected_labels) 
      
      TOP_X <- as.numeric(input$top_x)
      if (is.na(TOP_X) || TOP_X < 1) {
        showNotification("'Number of Top Features' must be a positive integer", type = "error")
        return()
      }
      
      if (TOP_X > nrow(data) ) {
        showNotification(paste0("'Number of Top Features' must be less than or equal to ", nrow(data)), type = "error")
        return()
      }
      
      # Generate the heatmap
      result <- plot_heatmap(data_subset, data, seq_subset, TOP_X, savedDatasetNameHeatmap(),
                             clustering_distance_rows, clustering_method_rows, 
                             show_column_names, show_row_names, cluster_rows,
                             show_row_dend, labels, enable_groups, groups, islog)
      
      heatmap_plot <- result$heatmap
      top_stats <- result$top_stats

      # Render the heatmap
      #output$heatmap_plot <- renderPlot({
      #  if (!is.null(heatmap_plot)) {
      #    draw(heatmap_plot)
      #  }
      #}, height = heatmap_height)

      pdf(file = NULL)  # Open a null PDF device to suppress output

      ht1 <- draw(heatmap_plot)
      makeInteractiveComplexHeatmap(input, output, session, ht1, "heatmap_interactive")
      
      # Render the table of top features
      output$heatmap_table <- DT::renderDataTable({
        DT::datatable(top_stats, options = list(pageLength = 20))
      })
      
      
      message(sample(quotes, 1))
      
    }
  })
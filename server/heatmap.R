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
          "Select Groups (in display order):", 
          choices = sort(unique(seq$group)),  # Still show sorted alphabetically
          selected = sort(unique(seq$group))[1],  # Default
          multiple = TRUE,
          width = "100%",
          selectize = TRUE  # This makes it searchable and orderable
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
      
      # # Extract column names from the selected dataset
      data_colnames <- colnames(data) # Substitude with data_sub
      columns <- c("heatmap_labels")
      for (column in columns) {
        # Update the 'identifier_column' select input with the new choices
        updateSelectInput(session, column, choices = data_colnames)
      }
    }
  })
  # Observer to update pickerInput choices based on selected grouping column
  observe({
    req(input$enable_grouping_heatmap, input$group_column_heatmap)
    
    # Get the current data
    if (input$select_heatmap_data == "Unsaved data") {
      data <- rv$tmpData
      req(!is.null(data))  # Make sure data exists
    } else {
      sd <- which(rv$choices %in% input$select_heatmap_data)
      req(length(sd) > 0)  # Make sure we found the dataset
      data <- rv$data[[sd]]
    }
    
    # Check if the selected grouping column exists in the data
    req(input$group_column_heatmap %in% colnames(data))
    
    # Get all unique values from the selected column
    all_values <- unique(data[[input$group_column_heatmap]])
    all_values <- sort(as.character(all_values[!is.na(all_values)]))
    
    # Debug message
    message(paste0("Updating pickerInput with ", length(all_values), " values from column: ", input$group_column_heatmap))
    
    # Update the pickerInput choices
    updatePickerInput(
      session = session,
      inputId = "selected_group_values",
      choices = all_values,
      selected = all_values  # Default: all selected
    )
  })
  
  # update when dataset changes
  observeEvent(input$select_heatmap_data, {
    req(input$enable_grouping_heatmap, input$group_column_heatmap)
    
    # Get the current data
    if (input$select_heatmap_data == "Unsaved data") {
      data <- rv$tmpData
      req(!is.null(data))
    } else {
      sd <- which(rv$choices %in% input$select_heatmap_data)
      req(length(sd) > 0)
      data <- rv$data[[sd]]
    }
    
    # Check if the selected grouping column exists
    if (input$group_column_heatmap %in% colnames(data)) {
      # Get all unique values from the selected column
      all_values <- unique(data[[input$group_column_heatmap]])
      all_values <- sort(as.character(all_values[!is.na(all_values)]))
      
      # Update the pickerInput choices
      updatePickerInput(
        session = session,
        inputId = "selected_group_values",
        choices = all_values,
        selected = all_values  # Default: all selected
      )
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
      
      # Check if "Name" column exists and is not empty
      if ("Name" %in% colnames(data)) {
        if (any(is.na(data[, "Name"]) | data[, "Name"] == "")) {
          sendSweetAlert(session, "Error",
                         "No names in Name column. Make sure features have names before generating heatmap.",
                         type = "error")
          return()
        }
      } else if ("name" %in% colnames(data)) {
        if (any(is.na(data[, "name"]) | data[, "name"] == "")) {
          sendSweetAlert(session, "Error",
                         "No names in name column. Make sure features have names before generating heatmap.",
                         type = "error")
          return()
        }
      } else {
        # If neither "Name" nor "name" column exists
        sendSweetAlert(session, "Error",
                       "No Name column found. Make sure features have names before generating heatmap.",
                       type = "error")
        return()
      }
      # Handle group selection if enabled
      if (input$select_groups_heatmap) {
        # User wants to select specific groups
        if (is.null(input$selected_groups_heatmap) || length(input$selected_groups_heatmap) == 0) {
          sendSweetAlert(session, "Error",
                         "Please select at least one group when 'Select Specific Groups' is enabled.",
                         type = "error")
          return()
        }
        
        # Get the selected groups (in the order they were selected)
        selected_groups <- input$selected_groups_heatmap
        
        # Filter seq_subset to only include selected groups
        seq_subset <- seq_subset[seq_subset$group %in% selected_groups, ]
        
        # SAFETY CHECK
        if (nrow(seq_subset) == 0) {
          sendSweetAlert(session, "Error",
                         "No samples match the selected groups.",
                         type = "error")
          return()
        }
        
        # Order the groups according to selection order
        seq_subset$group <- factor(seq_subset$group, levels = selected_groups)
        
      } else {
        # No specific group selection - include ALL groups
        # Just ensure groups are factors for consistent ordering
        all_groups <- unique(seq_subset$group)
        seq_subset$group <- factor(seq_subset$group, levels = sort(all_groups))
      }
      
      # Reorder seq_subset rows based on the factor order
      seq_subset <- seq_subset[order(seq_subset$group), ]
      
      # Filter data_subset columns based on the ordered seq_subset rownames
      data_subset <- data[, rownames(seq_subset), drop = FALSE]
      
      # SAFETY CHECK
      if (ncol(data_subset) == 0) {
        sendSweetAlert(session, "Error",
                       "No data columns remain after filtering.",
                       type = "error")
        return()
      }
      
      # SAFETY CHECK 2
      if (ncol(data_subset) == 0) {
        sendSweetAlert(session, "Error",
                       "No data columns remain after sample group filtering.",
                       type = "error")
        return()
      }
      #Filter features by selected group values if grouping is enabled
      if (input$enable_grouping_heatmap) {
        # Check if a grouping column is selected
        req(input$group_column_heatmap)
        
        # Check if any values are selected in the pickerInput
        if (is.null(input$selected_group_values) || length(input$selected_group_values) == 0) {
          sendSweetAlert(session, "Error",
                         "Please select at least one group value to display in the heatmap.",
                         type = "error")
          return()
        }
        
        groups <- input$group_column_heatmap
        
        # Check if the grouping column exists in the data
        if (!groups %in% colnames(data)) {
          showNotification(paste("Grouping column", groups, "not found in data."), type = "error")
          return()
        }
        
        # Get indices of rows (features) that match the selected group values
        feature_indices <- which(data[[groups]] %in% input$selected_group_values)
        
        if (length(feature_indices) == 0) {
          sendSweetAlert(session, "Error",
                         "No features match the selected group values.",
                         type = "error")
          return()
        }
        
        # Store original row count for messaging
        original_row_count <- nrow(data)
        
        # Filter data and data_subset to keep only selected features
        data <- data[feature_indices, , drop = FALSE]
        data_subset <- data_subset[feature_indices, , drop = FALSE]
        
        # SAFETY CHECK 3
        if (nrow(data) == 0) {
          sendSweetAlert(session, "Error",
                         "No features remain after group value filtering.",
                         type = "error")
          return()
        }
        
        # SAFETY CHECK 4
        if (nrow(data_subset) == 0) {
          sendSweetAlert(session, "Error",
                         "No data remains for heatmap after filtering.",
                         type = "error")
          return()
        }
        
        # Show success message about filtering
        showNotification(
          paste("Filtered to", length(feature_indices), "features based on group selection"),
          type = "message",
          duration = 3
        )
        
        message(paste0("Filtered to ", length(feature_indices), " features based on group selection"))
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
      
      # SAFETY CHECK 5
      if (length(selected_labels) == 0) {
        sendSweetAlert(session, "Error",
                       "No labels available for setting row names.",
                       type = "error")
        return()
      }
      
      if (nrow(data_subset) == 0) {
        sendSweetAlert(session, "Error",
                       "data_subset is empty before setting row names.",
                       type = "error")
        return()
      }
      
      rownames(data_subset) <- make.unique(selected_labels)
      rownames(data) <- make.unique(selected_labels) 
      
      # NEW: Update row names pickerInput with the FINAL row names
      if (input$show_row_names) {
        final_row_names <- rownames(data)
        updatePickerInput(
          session = session,
          inputId = "selected_row_names",
          choices = final_row_names,
          selected = final_row_names  # Default to all selected
        )
      }
      
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
      
      # NEW: Check if rows were removed due to missing values
      if (result$rows_removed > 0) {
        if (nrow(top_stats) < TOP_X) {
          showNotification(
            paste("Showing", nrow(top_stats), "of", TOP_X, "requested features -", 
                  result$rows_removed, "features removed due to missing values"),
            type = "warning",
            duration = 5
          )
        } else {
          showNotification(
            paste("Removed", result$rows_removed, "features with missing values"),
            type = "message",
            duration = 3
          )
        }
      }
      
      # NEW: Filter the top_stats table based on selected row names (if any)
      if (input$show_row_names && !is.null(input$selected_row_names) && length(input$selected_row_names) > 0) {
        top_stats <- top_stats[top_stats$Feature %in% input$selected_row_names, , drop = FALSE]
        
        if (nrow(top_stats) == 0) {
          showNotification("No data left after filtering by selected row names.", type = "warning")
        } else {
          showNotification(
            paste("Showing", nrow(top_stats), "rows based on selection"),
            type = "message",
            duration = 2
          )
        }
      }

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
  ########## NEW OBSERVER FOR top_x CHANGES ##########
  observeEvent(input$top_x, {
    # Only trigger if a heatmap has already been generated
    if (!is.null(rv$last_heatmap)) {
      shinyjs::click("run_heatmap")
    }
  })
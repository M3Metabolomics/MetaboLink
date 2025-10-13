  # Generate Volcano Plot
  
  observeEvent(input$select_volcano_data, {
    req(input$select_volcano_data) # Ensure a dataset is selected
    
    if (!is.null(rv$activeFile)) {
      if (input$select_volcano_data == "Unsaved data") {
        data <- rv$tmpData  # Use the temporary data
      } else {
        # Get the index of the selected dataset
        sd <- which(rv$choices %in% input$select_volcano_data)
        data <- rv$data[[sd]]  # Retrieve the selected dataset
      }
      
      # Extract column names from the selected dataset
      data_colnames <- colnames(data)
      
      columns <- c("volcano_labels")
      
      for (column in columns) {
        # Update the 'identifier_column' select input with the new choices
        updateSelectInput(session, column, choices = data_colnames)
      }
    }
  })
  # Define the reactive value at the top of the server 
  savedDatasetNameVolcano <- reactiveVal("Volcano Plot: ")
  # Observe the input for the volcano title and update the reactive value
  observe({
    savedDatasetNameVolcano(input$volcano_title)
    
    output$displayName <- renderText({
      paste("Current Heatmap Title:", savedDatasetNameVolcano())
    })
  })
  
  output$feature_selection_ui_volcano <- renderUI({
    if (!input$enable_feature_selection) return(NULL)  # Only show if checkbox is checked
    req(input$select_volcano_data)  # Ensure a dataset is selected
    
    if (!is.null(rv$activeFile)) {
      if (input$select_volcano_data == "Unsaved data") {
        data <- rv$tmpData  # Use the temporary data
        seq <- rv$tmpSequence  # Use the temporary sequence
      } else {
        # Get the index of the selected dataset
        sd <- which(rv$choices %in% input$select_volcano_data)
        data <- rv$data[[sd]]  # Retrieve the selected dataset
        seq <- rv$sequence[[sd]]  # Retrieve the selected sequence
      }
      
      # Filter sequence to include only "Sample" rows
      seq <- seq[!seq[, "labels"] %in% c("Sample", "QC"), ]
      data <- data[, rownames(seq), drop = FALSE]  # Filter data based on seq row names
      
      columns <- colnames(data)  # Extract column names
      
      default_val <- if ("refmet_name" %in% columns) {
        "refmet_name"
      } else {
        columns[1]  # Default to first column if "refmet_name" is not found
      }
      
      # **Fix: Wrap in `tagList()` so both UI elements render correctly**
      tagList(
        selectInput(
          "volcano_feature_column",  # Input ID for selecting feature column
          "Select Feature Column:",
          choices = columns,
          selected = default_val,
          width = "100%"
        ),
        
        selectizeInput(
          "selected_features_volcano",
          "Select Features:",
          choices = NULL,  # Choices will be updated dynamically
          multiple = TRUE,
          options = list(
            placeholder = "Search & Select Features",
            maxOptions = 100  # Show only 100 at a time
          )
        )
      )
    }
  })
  observeEvent(input$volcano_feature_column, {
    req(input$volcano_feature_column, input$select_volcano_data)  # Ensure valid inputs
    
    if (!is.null(rv$activeFile)) {
      if (input$select_volcano_data == "Unsaved data") {
        data <- rv$tmpData
      } else {
        sd <- which(rv$choices %in% input$select_volcano_data)
        data <- rv$data[[sd]]
      }
      
      # Get the unique feature names from the selected column
      feature_choices <- unique(data[[input$volcano_feature_column]])
      
      updateSelectizeInput(
        session,
        "selected_features_volcano",
        choices = feature_choices,
        selected = NULL,  # Reset selection
        server = TRUE  # **Enable Server-Side Processing**
      )
    }
  })
  
  output$group_selection_ui_volcano <- renderUI({
    if (!input$enable_group_selection) return(NULL)
    req(input$select_volcano_data)
    
    # Retrieve data (this example is based on your existing code)
    if (!is.null(rv$activeFile)) {
      if (input$select_volcano_data == "Unsaved data") {
        data <- rv$tmpData  
        seq <- rv$tmpSequence  
      } else {
        sd <- which(rv$choices %in% input$select_volcano_data)
        data <- rv$data[[sd]]
        seq <- rv$sequence[[sd]]
      }
      
      seq <- seq[!seq[, "labels"] %in% c("Sample", "QC"), ]
      data <- data[, rownames(seq), drop = FALSE]
      columns <- colnames(data)
      
      default_val <- if ("sub_class" %in% columns) {
        "sub_class"
      } else if ("Lipid.Abbreviation" %in% columns) {
        "Lipid.Abbreviation"
      } else {
        columns[1]
      }
    }
    
    tagList(
      selectInput(
        "volcano_group_column",
        "Select Group Column:",
        choices = columns,
        selected = default_val,
        width = "100%"
      ),
      selectizeInput(
        "selected_group_volcano",
        "Select Group:",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "Search & Select Groups")
      )
    )
  })
  
  darken_color <- function(color, factor = 0.7) {
    rgb_val <- grDevices::col2rgb(color)
    dark_rgb <- pmax(rgb_val * factor, 0)
    rgb(t(dark_rgb), maxColorValue = 255)
  }
  
  output$group_color_ui <- renderUI({
    # Only show if group selection is enabled and at least one group is selected
    if (!input$enable_group_selection) return(NULL)
    req(input$selected_group_volcano)
    
    selected_groups <- input$selected_group_volcano
    n_groups <- length(selected_groups)
    default_colors <- hue_pal()(n_groups)

    # For each selected group, create two colourInput widgets:
    ui_list <- lapply(seq_along(selected_groups), function(i) {
      grp <- selected_groups[i]
      grp_id <- make.names(grp)
      
      fill_default <- default_colors[i]
      outline_default <- darken_color(fill_default)
      
      tagList(
        h4(paste("Group:", grp)),
        fluidRow(
          column(
            width = 6,
            colourInput(
              inputId = paste0("color_", grp_id, "_fill"),
              label = "Fill:",
              value = fill_default  # Default fill color
            )
          ),
          column(
            width = 6,
            colourInput(
              inputId = paste0("color_", grp_id, "_outline"),
              label = "Outline:",
              value = outline_default  # Default outline color
            )
          )
        )
      )
    })
    
    tagList(ui_list)
  })
  observeEvent(input$volcano_group_column, {
    req(input$volcano_group_column,
        input$select_volcano_data)  # Ensure valid inputs
    
    if (!is.null(rv$activeFile)) {
      if (input$select_volcano_data == "Unsaved data") {
        data <- rv$tmpData
      } else {
        sd <- which(rv$choices %in% input$select_volcano_data)
        data <- rv$data[[sd]]
      }
      
      # Ensure the selected column exists in data before accessing it
      if (!input$volcano_group_column %in% colnames(data)) return()
      
      unique_groups <- unique(data[[input$volcano_group_column]])
      
      updateSelectizeInput(
        session,
        "selected_group_volcano",
        choices = unique_groups,
        selected = unique_groups[2],  # Reset selection
        server = TRUE  # Enable server-side processing for large lists
      )
      
    }
  })
  
  output$parameter_selection_ui_volcano <- renderUI({
    if (isTRUE(input$select_parameter_volcano)) {
      fluidRow(
        column(
          width = 6,
          numericInput("x_param", "X axis limit", value = 5)
        ),
        column(
          width = 6,
          numericInput("y_param", "Y axis limit", value = 5)
        )
      )
    }
  })
  
  observeEvent(input$run_volcano_plot, {
    req(
      input$select_volcano_data,
      input$volcano_labels,
      input$group1_vol, 
      input$group2_vol,
      input$log2fc_threshold, 
      input$pval_threshold,
      input$color_up_fill,
      input$color_up_outline,
      input$color_down_fill,
      input$color_down_outline,
      input$color_ns_fill,
      input$color_ns_outline
    )

    if (!is.null(rv$activeFile)) {
      if (input$select_volcano_data == "Unsaved data") {
        data <- rv$tmpData  # Use the temporary data
        seq <- rv$tmpSequence  # Use the temporary sequence
        dataset_name <- "Unsaved data"
      } else {
        # Get the index of the selected dataset
        sd <- which(rv$choices %in% input$select_volcano_data)
        data <- rv$data[[sd]]  # Retrieve the selected dataset
        seq <- rv$sequence[[sd]]  # Retrieve the selected sequence
        dataset_name <- names(rv$data)[sd]  # Retrieve dataset name
      }
      
      label_column <- input$volcano_labels
      numerator <- input$group1_vol
      denominator <- input$group2_vol
      
      # make an error check that numerator and denominator are not the same
      if (numerator == denominator) {
        showNotification(
          "The numerator and denominator groups must be different!",
          type = "warning"
        )
        return(NULL) # return from the function/observe and don't proceed
      }
      
      log2FC_tresh <- log2(input$log2fc_threshold)
      pAdjustMethod <- input$pAdjustMethod_volcano
      pval_tresh <- input$pval_threshold
      pval_col <- input$pval_col_volcano
      
      fill_up <- input$color_up_fill
      outline_up <- input$color_up_outline
      fill_down <- input$color_down_fill
      outline_down <- input$color_down_outline
      fill_ns <- input$color_ns_fill
      outline_ns <- input$color_ns_outline
      
      show_legend <- input$show_legend_volcano
      
      x_param <- input$x_param
      y_param <- input$y_param
      apply_axis_limits <- input$select_parameter_volcano
      
      message(paste0("dataset_name: ", dataset_name))
      message(paste0("label column: ", label_column))
      message(paste0("numerator: ", numerator))
      message(paste0("denominator: ", denominator))
      message(paste0("log2FC input: ", 2^(log2FC_tresh)))
      message(paste0("log2FC tresh: ", log2FC_tresh))
      message(paste0("pval tresh: ", pval_tresh))
      message(paste0("show legend: ", show_legend))
      
      message(paste0("Axis limits: ", apply_axis_limits))
      
      if (input$select_parameter_volcano) {
        message(paste0("X axis limit: ", x_param))
        message(paste0("Y axis limit: ", y_param))
      }
    
      enable_feature_selection <- input$enable_feature_selection
      message(paste0("Feature Selection Enabled: ", enable_feature_selection))
      if (enable_feature_selection) {
        available_features <- input$selected_features_volcano
        message("Available Features:")
        print(available_features)
      }
      
      enable_group_selection <- input$enable_group_selection
      message(paste0("Group Selection Enabled: ", enable_group_selection))
      if (enable_group_selection) {
        available_groups <- input$selected_group_volcano
        message("Available Groups:")
        print(available_groups)
      }
      
      
      if (input$enable_group_selection && !is.null(input$selected_group_volcano)) {
        selected_groups <- input$selected_group_volcano
        group_color_df <- do.call(rbind, lapply(selected_groups, function(grp) {
          grp_id <- make.names(grp)
          data.frame(
            Group   = grp,
            Fill    = input[[paste0("color_", grp_id, "_fill")]],
            Outline = input[[paste0("color_", grp_id, "_outline")]],
            stringsAsFactors = FALSE
          )
        }))
        
        # Mapping each group to its fill and outline colors.
        print(group_color_df)
      }
      
      if (enable_feature_selection) {
        if (length(available_features) == 0) {
          showNotification(
            "Select at least one feature for the volcano plot.",
            type = "warning"
          )
          return(NULL)  # Stop execution if no features are selected
        }
        
        showNotification("This feature is of limited use and may not work as expected. ",
                 "As of 5/3-2025 this feature is still under development.", type = "message")
        
        # Check if feature column exists
        if (input$volcano_feature_column %in% colnames(data)) {
          # Filter data to include only selected features
          feature_sub <- data[data[[input$volcano_feature_column]] %in% available_features, , drop = FALSE]
          
          # Print only the first few rows of the selected feature subset
          print("Subset of Selected Features:")
          print(head(feature_sub,2))
        } else {
          showNotification("Feature column not found in dataset.", type = "error")
        }
      }
      
      if (enable_group_selection) {
        if (length(available_groups) == 0) {
          showNotification(
            "Select at least one group for the volcano plot.",
            type = "warning"
          )
          return(NULL)  # Stop execution if no groups are selected
        }
        
        # Check if group column exists
        if (input$volcano_group_column %in% colnames(data)) {
          # Filter data to include only selected groups
          group_sub <- data[data[[input$volcano_group_column]] %in% available_groups, , drop = FALSE]
          
          # Print only the first few rows of the selected group subset
          print("Subset of Selected Groups:")
          print(head(group_sub,2))
        } else {
          showNotification("Group column not found in dataset.", type = "error")
        }
      }
      
      seq_subset <- seq[seq[, "labels"] %in% c("Sample"), ]  # Restrict to "Sample" rows
      data_subset <- data[, c(rownames(seq_subset)), drop = FALSE]  # Use row names of seq_subset to filter columns
      
      selected_labels <- as.character(data[[label_column]])
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
      
      stat_results <- calculate_stats(data_subset, seq_subset, adjust.method = pAdjustMethod)
      
      target_contrast   <- paste0(numerator, "_vs_", denominator)
      reversed_contrast <- paste0(denominator, "_vs_", numerator)
      
      # If  contrast is present in stat_results:
      if (target_contrast %in% stat_results$Contrast) {
        # Just subset
        sub_df <- subset(stat_results, Contrast == target_contrast)
      } else if (reversed_contrast %in% stat_results$Contrast) {
        # Subset, then flip the log2FC & FC
        sub_df <- subset(stat_results, Contrast == reversed_contrast)
        sub_df$log2FC <- -sub_df$log2FC
        sub_df$FC     <- 1 / sub_df$FC
        
        # Rename the contrast column to reflect the new direction
        sub_df$Contrast <- target_contrast
      } else {
        # No matching contrast found; handle how you like (warn user, or return empty)
        warning(
          paste0("No matching contrast found for '", numerator, " vs ", denominator, "'. ",
                 "Available contrasts are: ", paste(unique(stat_results$Contrast), collapse=", "))
        )
        sub_df <- data.frame()
      }
      
      data <- data %>%
        rownames_to_column("Feature_ID")
      
      sub_df <- sub_df %>%
        rownames_to_column("Feature_ID") %>%
        relocate(Feature_ID, .before = "Contrast")
      
      # make the rownames as the rownames of the sub_df
      rownames(sub_df) <- sub_df$Feature_ID
      sub_df$Feature_ID <- gsub("^[^.]+\\.", "", sub_df$Feature_ID)
      
      if (input$enable_group_selection) {
        sub_df <- sub_df %>%
          left_join(
            data %>% select(Feature_ID, !!sym(input$volcano_group_column)),
            by = "Feature_ID"
          ) %>%
          rename(Group = !!sym(input$volcano_group_column)) %>%
          relocate(Group, .before = "Contrast")
      }
      
      # assign the sub_df a name for debugging
      volcano_df_name <- paste0(dataset_name, "_volcano_df")
      assign(volcano_df_name, sub_df)
      
      pvp <- reactive(volcano_plot(
        sub_df,savedDatasetNameVolcano(),
        log2FC_tresh, pval_tresh,
        fill_up, outline_up,
        fill_down, outline_down,
        fill_ns, outline_ns,
        enable_feature_selection, enable_group_selection,
        available_features, available_groups,
        group_color_df,x_param,y_param, apply_axis_limits,
        pval_col))
      
      # Render the volcano plot
      output$volcano_plot <- renderPlotly({
        pvp()$plot
      })
      
      # Render the table of top features
      output$volcano_table <- DT::renderDataTable({
        df <- pvp()$df
        df <- df %>%
          select(-c(hover_text,Contrast))
        DT::datatable(df, 
                      rownames = FALSE,
                      options = list(scrollX = TRUE,
                                     pageLength = 20,
                                     autoWidth = TRUE))
      })
      
      message(sample(quotes, 1))
    }
  })
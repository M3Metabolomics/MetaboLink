#################
# Lipid Heatmap #
#################
# Values is used for message display before and after data process
values <- reactiveValues(
  runProcessClicked = FALSE,
  removed_data = NULL,
  grouped_data_frames = NULL,
  lipid_names = NULL,
  cleaning_stats = NULL  # <--- NEW: Stores the counts for the summary
)


# Create the inputs once; we will update their choices/selected later
output$select_group_ui_heatmap <- renderUI({
  req(isTRUE(values$runProcessClicked))   # <- hides UI until button is clicked
  
  column(
    width = 12,
    tagList(
      selectInput("selected_group_for_numerator",  "Select Group for numerator:",  choices = character(0)),
      selectInput("selected_group_for_denominator","Select Group for denominator:", choices = character(0))
    )
  )
})

# Remember last chosen groups across reprocessing / dataset switches
values$last_num <- NULL
values$last_den <- NULL

observeEvent(input$selected_group_for_numerator,  ignoreInit = TRUE, {
  values$last_num <- input$selected_group_for_numerator
})
observeEvent(input$selected_group_for_denominator, ignoreInit = TRUE, {
  values$last_den <- input$selected_group_for_denominator
})






# Renders the selector (hidden until a dataset is available)
output$heatmap_name_col_ui <- renderUI({
  req(rv$activeFile, rv$data[[rv$activeFile]])
  this_data <- rv$data[[rv$activeFile]]
  cols <- colnames(this_data)
  
  
  # choose a default: prefer a name-ish column if present
  default_val <- if ("Compound_Name" %in% cols) {
    "Compound_Name"
  } else if ("Name" %in% cols) {
    "Name"
  } else if ("Original annotation" %in% cols) {
    "Original annotation"
  } else {
    cols[1]
  }
  
  tagList(
    br(),
    selectInput(
      "heatmap_name_col",
      "Select column with lipid names:",
      choices  = cols,
      selected = default_val,
      width = "100%"
    )
  )
})

# If  active dataset can change, keep the choices in sync:
observeEvent(rv$activeFile, {
  req(rv$activeFile, rv$data[[rv$activeFile]])
  cols <- colnames(rv$data[[rv$activeFile]])
  
  # try to keep the previous selection if it still exists
  sel <- if (!is.null(input$heatmap_name_col) && input$heatmap_name_col %in% cols) {
    input$heatmap_name_col
  } else if ("Compound_Name" %in% cols) {
    "Compound_Name"
  } else if ("Name" %in% cols) {
    "Name"
  } else if ("Original annotation" %in% cols) {
    "Original annotation"
  } else {
    cols[1]
  }
  
  updateSelectInput(session, "heatmap_name_col", choices = cols, selected = sel)
})


# When button clicked in interface, all the following will be processed
observeEvent(input$run_process, {
  
  
  # --- FIX START: Check if Sequence Exists ---
  # We check three things to be safe:
  # 1. Is there an active file selected? (rv$activeFile)
  # 2. Is the sequence list long enough to contain this file?
  # 3. Is the sequence data for this file actually there (not NULL)?
  
  has_sequence <- FALSE
  if (!is.null(rv$activeFile) && length(rv$sequence) >= rv$activeFile) {
    if (!is.null(rv$sequence[[rv$activeFile]])) {
      has_sequence <- TRUE
    }
  }
  
  if (!has_sequence) {
    # Option A: Use shinyalert (Prettier, requires useShinyalert() in UI)
    shinyalert::shinyalert(
      title = "Missing Sequence File",
      text = "Please upload a sequence file (metafile) before running data processing.",
      type = "warning"
    )
    
    
    
    return() # <--- CRITICAL: This stops the code here so it doesn't crash below!
  }
  # --- FIX END ---
  
  
  
  values$runProcessClicked <- TRUE
  
  sequence <- rv$sequence[[rv$activeFile]]
  data     <- rv$data[[rv$activeFile]]
  
  # --- STEP 1: Handle the User's Column Selection ---
  req(input$heatmap_name_col)
  name_col <- input$heatmap_name_col
  
  # Relocate and Rename
  if (name_col %in% colnames(data)) {
    data <- dplyr::relocate(data, dplyr::all_of(name_col), .before = 1)
  }
  colnames(data)[1] <- "Compound_Name"
  
  # 1. CAPTURE START COUNT
  n_start <- nrow(data)
  
  ###############
  # Data cleaning
  ###############
  
  # Apply formatting
  data[, 1] <- sapply(data[, 1], format_strings)
  
  # --- FIX: Normalize Lipid Class Case (pc -> PC) ---
  # This finds the text before the '(' and forces it to UPPERCASE.
  data[, 1] <- sub("^([^(]+)", "\\U\\1", data[, 1], perl = TRUE)
  
  data[, 1] <- sapply(data[, 1], make_total_name, only_if_multiple = FALSE)
  
  # 2. CAPTURE & REMOVE INVALID PATTERNS
  # Store the rows that WILL be removed (for the table)
  values$removed_data <- remove_patterned_rows(data) 
  n_invalid_format <- nrow(values$removed_data)
  
  # Filter rows (Keep only valid X(C:D) format)
  data <- filter_data_by_pattern(data)
  n_post_clean <- nrow(data)
  
  # Initialize info object
  merged_data_info <- NULL
  
  # --- STEP 2: Handle Dataset Selection (Original vs Merged) ---
  n_merged_loss <- 0 # Default to 0
  
  if (input$selected_dataset == "original") {
    data <- unique_compound_names(data)
  } else if (input$selected_dataset == "merged") {
    merged_data_info <- merged_info_function(data)
    
    # Identify sample columns
    sample_cols_indices <- which(sequence[, 'labels'] == "Sample")
    data <- data[, c(1, sample_cols_indices), drop = FALSE]
    sequence <- sequence[sequence[, 'labels'] %in% c("Name", "Sample"), ]
    
    # Run Merge
    n_before_merge <- nrow(data)
    data <- merge_duplicates(data)
    n_after_merge <- nrow(data)
    
    # Calculate how many rows were collapsed
    n_merged_loss <- n_before_merge - n_after_merge
  }
  
  # 3. CAPTURE FINAL COUNT
  n_final <- nrow(data)
  
  # SAVE STATS TO VALUES (So we can show them in the modal)
  values$cleaning_stats <- list(
    start = n_start,
    invalid = n_invalid_format,
    valid = n_post_clean,
    merged_loss = n_merged_loss,
    final = n_final
  )
  
  
  # The following is used in the tab: 'Lipid summary'.
  #grouped_samples <- process_lipid_data(sequence, data) # not used?
  output$groups_table <- renderTable({
    
    # Extract the 'Sample' labels and corresponding 'group' from 'sequence'
    sample_rows <- sequence[sequence$labels == "Sample", ]
    unique_groups <- unique(sample_rows$group)
    
    # Create the dataframe to be displayed as a table
    lipid_df_processed_data <- data.frame(
      Group = unique_groups,
      Samples = sapply(unique_groups, function(group) {
        sample_identifiers <- rownames(sample_rows)[sample_rows$group == group]
        paste(sample_identifiers, collapse = ", ")
      })
    )
    # Return the data frame to be rendered as a table
    lipid_df_processed_data
  })
  
  
  
  # Heatmap input selection  
  observeEvent(input$run_process, {
    processed_results <- process_lipid_data(sequence, data)
    values$grouped_data_frames <- create_grouped_data_frames(sequence, data)
    
    # --- FIX START: Stop if no matching data found ---
    if (is.null(values$grouped_data_frames) || length(values$grouped_data_frames) == 0) {
      shinyalert::shinyalert(
        title = "Data Mismatch",
        text = "Could not group data. Please check if you have uploade Sequences file, or if your Sequence file matches the column names in your Data file.",
        type = "error"
      )
      values$runProcessClicked <- FALSE # Reset the button state
      return() # Stop here to prevent the crash
    }
    # --- FIX END ---
    
    compound_names <- data[[1]]  # Extract the first column which contains compound names
    
    # Assuming that each grouped data frame has rows in the same order as "data"
    for (i in seq_along(values$grouped_data_frames)) {
      values$grouped_data_frames[[i]] <- cbind(Compound_Name = compound_names, values$grouped_data_frames[[i]])
    }
    
    # Extract unique group names from sequence
    unique_group_names <- unique(sequence[sequence$labels == "Sample", "group"])
    
    # Check if lengths of group names and grouped data frames match
    if (length(unique_group_names) == length(values$grouped_data_frames)) {
      # Apply the actual group names from the sequence file
      names(values$grouped_data_frames) <- unique_group_names
    } else {
      # If there's a mismatch, fallback to naming with numbers as before
      names(values$grouped_data_frames) <- paste("Group", seq_along(values$grouped_data_frames))
    }
    
    choices <- names(values$grouped_data_frames)
    
    # fallbacks if previous selection isn’t available in this dataset
    sel_num <- if (!is.null(values$last_num) && values$last_num %in% choices) {
      values$last_num
    } else {
      choices[1]
    }
    
    sel_den <- if (!is.null(values$last_den) && values$last_den %in% choices) {
      values$last_den
    } else {
      choices[ pmin(2, length(choices)) ]  # second choice if it exists, else first
    }
    
    updateSelectInput(session, "selected_group_for_numerator",
                      choices = choices, selected = sel_num)
    updateSelectInput(session, "selected_group_for_denominator",
                      choices = choices, selected = sel_den)
    
    
    
    
    
    
    # Create interactive table for selected numerator group. Displayed at the starting page of the 'Lipid Heatmap'
    output$numerator_group_table <- DT::renderDataTable({
      req(input$selected_group_for_numerator, values$grouped_data_frames) 
      # Create a copy of the data for display purposes
      display_data <- values$grouped_data_frames[[input$selected_group_for_numerator]]
      
      DT::datatable(
        display_data,  
        options = list(scrollX = TRUE)  # Enable horizontal scrolling
      )
    })
    
    # Create interactive table for selected denominator group
    output$denominator_group_table <- DT::renderDataTable({
      req(input$selected_group_for_denominator, values$grouped_data_frames)  
      display_data <- values$grouped_data_frames[[input$selected_group_for_denominator]]
      
      DT::datatable(
        display_data,  
        options = list(scrollX = TRUE)  
      )
    })
    
    
    output$select_lipid_ui <- renderUI({
      # Extract lipid classes from the first column of the current data
      values$lipid_names <- group_lipids_by_group(data)
      
      # Available choices
      choices <- c("All", sort(unique(values$lipid_names$group)))
      
      # Try to restore previous selection, but only keep items that still exist
      selected <- values$last_lipids
      if (!is.null(selected)) {
        selected <- intersect(selected, choices)
        # fallback if nothing intersects
        if (length(selected) == 0) selected <- "All"
      } else {
        selected <- "All"
      }
      
      selectizeInput(
        "selected_lipid", "Select lipid(s) to display:",
        choices  = choices,
        selected = selected,
        multiple = TRUE,
        options  = list(
          placeholder = 'Choose lipids...',
          plugins = list('remove_button')  # nice UX; optional
        )
      )
    })
    
    
    # Reactive expression to track the number of selected lipids or the "All" selection
    selected_lipid_count <- reactive({
      # If "All" is selected, we could set this to a value that causes the default text size to be used
      if ("All" %in% input$selected_lipid) {
        return(Inf)  # 'Inf' is used here as a flag for "All"
      } else {
        return(length(input$selected_lipid))
      }
    })
    
    
    
    reactiveP_value <- reactive({
      req(input$selected_group_for_numerator, input$selected_group_for_denominator)
      
      # FIX: Use values$grouped_data_frames
      numerator_data <- values$grouped_data_frames[[input$selected_group_for_numerator]]
      denominator_data <- values$grouped_data_frames[[input$selected_group_for_denominator]]
      
      # Ensure there is data to work with
      if (nrow(numerator_data) == 0 || nrow(denominator_data) == 0) {
        return(NULL)
      }
      
      # Initialize a vector to store the p-values
      p_values <- numeric(nrow(numerator_data))
      
      # Loop through each row to perform the t-test
      for (i in 1:nrow(numerator_data)) {
        # Extract the numerical values for numerator and denominator, excluding the first column
        num_values <- numerator_data[i, -1]
        denom_values <- denominator_data[i, -1]
        
        # Check if data is constant or contains NA values
        if (length(unique(num_values)) == 1 || length(unique(denom_values)) == 1 ||
            any(is.na(num_values)) || any(is.na(denom_values))) {
          p_values[i] <- NA  # Assign NA or another appropriate value
        } else {
          # Perform the t-test
          # FIX: Wrapped in try() to prevent crashing on single errors
          try({
            t_test_result <- t.test(num_values, denom_values)
            p_values[i] <- t_test_result$p.value
          }, silent = TRUE)
        }
      }
      
      # Apply Benjamini-Hochberg correction to p-values
      adjusted_p_values <- p.adjust(p_values, method = "BH")
      
      # Create a new data frame with 'Compound_Name', 'p_value', and 'padj'
      p_value_data <- data.frame(
        Compound_Name = numerator_data$Compound_Name,  
        p_value = p_values,
        padj = adjusted_p_values
      )
      
      # Filter p-values on adjusted_p_values based on ui input
      p_value_data <- p_value_data[!is.na(p_value_data$padj) & p_value_data$padj < input$p_value_adj, ]
      return(p_value_data)
    })
    
    # Reactive expression to calculate logFC
    reactiveLogFC <- reactive({
      # The required data input for the data handling. 
      req(input$selected_group_for_numerator, input$selected_group_for_denominator)
      
      # Define data input, makes it more readable 
      # FIX: Use values$grouped_data_frames
      numerator_data <- values$grouped_data_frames[[input$selected_group_for_numerator]]
      denominator_data <- values$grouped_data_frames[[input$selected_group_for_denominator]]
      
      # Removes the fir column of the data, as it is not needed for the calculation.
      numerator_data <- numerator_data[, -1]
      denominator_data <- denominator_data[, -1]
      
      numerator_data_means <- rowMeans(numerator_data, na.rm = TRUE) 
      denominator_data_means <- rowMeans(denominator_data, na.rm = TRUE) 
      numerator_data_means <- data.frame(numerator_data_means)
      denominator_data_means <- data.frame(denominator_data_means)
      
      # Extract the compound names, to add it to the data frame logFC, making sure they are sorted by compound names. 
      compound_names <- data[[1]]
      
      # Calculate logFC
      logFC_data <- log2((numerator_data_means + 1e-6) / (denominator_data_means + 1e-6))
      # Rename a single column
      colnames(logFC_data)[colnames(logFC_data) == "numerator_data_means"] <- "logFC"
      
      logFC <- data.frame(Compound_Name = compound_names, logFC = logFC_data)
      
      
      # Continue filtering based on lipid selection
      if (!"All" %in% input$selected_lipid) {
        # FIX: Use values$lipid_names
        logFC <- logFC[values$lipid_names$group %in% input$selected_lipid, ]
      }
      
      filtered_data <- logFC
      return(filtered_data) # Used for Heatmap display 
    })
    
    
    
    ##### Render UI for different thresholds within the app. 
    
    # Render UI for maximum p-value input
    output$p_value_max_ui <- renderUI({
      numericInput("p_value_max", 
                   "Maximum p-value:", 
                   value = 1, 
                   min = 0, 
                   step = 0.01)
    })
    
    # Render UI for maximum p-value input
    output$p_value_adj <- renderUI({
      numericInput("p_value_adj", 
                   "Maximum p-value_adj:", 
                   value = 1, 
                   min = 0, 
                   step = 0.01)
    })
    
    # Render UI for logFC input
    output$logFC_input_ui <- renderUI({
      tagList(
        numericInput("logFC_input", 
                     "Enter logFC threshold:", 
                     value = 0,
                     min = 0)
      )
    })
    
    # Add this in the UI section where other inputs are rendered
    output$min_lipids_per_class_ui <- renderUI({
      numericInput("min_lipids_per_class", 
                   "Minimum number of lipids per class:", 
                   value = 2, 
                   min = 1, 
                   step = 1)
    })
    
    
    
    #### Filtration within the app
    
    # Reactive expression to filter data based on p-value and logFC thresholds,
    # plus the amount of lipids within their class.
    #### Filtration within the app
    
    # Reactive expression to filter data based on:
    # - p-value threshold
    # - logFC threshold
    # - Carbon chain range
    # - Double bonds range
    # - minimum lipids per FacetClass
    reactiveFilteredData <- reactive({
      logFC_data   <- reactiveLogFC()
      p_value_data <- reactiveP_value()
      
      # Join on Compound_Name (padj filtering already applied inside reactiveP_value)
      filtered_data <- dplyr::inner_join(logFC_data, p_value_data, by = "Compound_Name")
      
      # Basic thresholds: p-value & logFC
      filtered_data <- filtered_data %>%
        dplyr::filter(!is.na(.data$p_value) & .data$p_value <= input$p_value_max) %>%
        dplyr::filter(!is.na(.data$logFC)   & abs(.data$logFC) >= input$logFC_input)
      
      # If nothing left -> return early
      if (nrow(filtered_data) == 0) {
        return(filtered_data)
      }
      
      # ---- Name mapping: Class / Carbon / Double / FacetClass ----
      nm <- compute_names_mapping(filtered_data$Compound_Name)
      
      # If mapping failed for some reason, just return what we have
      if (nrow(nm) == 0) {
        return(filtered_data[0, , drop = FALSE])
      }
      
      # FacetClass: e.g. "TG_1(58:8)" -> "TG_1"
      nm$FacetClass <- sub("\\s*\\(.*$", "", nm$Compound_Name)
      
      # make Carbon / Double numeric for filtering
      nm$Carbon <- suppressWarnings(as.integer(nm$Carbon))
      nm$Double <- suppressWarnings(as.integer(nm$Double))
      
      # ---- Carbon chain range filter (if slider exists) ----
      if (!is.null(input$carbon_range) && length(input$carbon_range) == 2) {
        nm <- nm[
          !is.na(nm$Carbon) &
            nm$Carbon >= input$carbon_range[1] &
            nm$Carbon <= input$carbon_range[2],
          ,
          drop = FALSE
        ]
      }
      
      # ---- Double bonds range filter (if slider exists) ----
      if (!is.null(input$double_range) && length(input$double_range) == 2) {
        nm <- nm[
          !is.na(nm$Double) &
            nm$Double >= input$double_range[1] &
            nm$Double <= input$double_range[2],
          ,
          drop = FALSE
        ]
      }
      
      # If all lipids were removed by Carbon/Double filtering
      if (nrow(nm) == 0) {
        return(filtered_data[0, , drop = FALSE])
      }
      
      # ---- Facet-wise min lipids per class (using FacetClass from nm) ----
      facet_counts_df <- nm %>%
        dplyr::filter(!is.na(.data$FacetClass)) %>%
        dplyr::count(.data$FacetClass, name = "Count")
      
      facets_to_keep <- facet_counts_df$FacetClass[
        facet_counts_df$Count >= input$min_lipids_per_class
      ]
      
      # Keep only compounds in FacetClasses that meet the threshold
      nm <- nm[
        !is.na(nm$FacetClass) & nm$FacetClass %in% facets_to_keep,
        ,
        drop = FALSE
      ]
      
      if (nrow(nm) == 0) {
        return(filtered_data[0, , drop = FALSE])
      }
      
      # ---- Finally, restrict filtered_data to the surviving Compound_Names ----
      filtered_data <- filtered_data[
        filtered_data$Compound_Name %in% nm$Compound_Name,
        ,
        drop = FALSE
      ]
      
      filtered_data
    })
    
    
    # p-values for ALL lipids (no filtering); used only for drawing outlines
    reactiveP_value_all <- reactive({
      req(input$selected_group_for_numerator, input$selected_group_for_denominator)
      
      # FIX: Use values$grouped_data_frames
      numerator_data   <- values$grouped_data_frames[[input$selected_group_for_numerator]]
      denominator_data <- values$grouped_data_frames[[input$selected_group_for_denominator]]
      
      if (nrow(numerator_data) == 0 || nrow(denominator_data) == 0) return(NULL)
      
      p_values <- numeric(nrow(numerator_data))
      
      for (i in 1:nrow(numerator_data)) {
        num_values   <- numerator_data[i, -1]
        denom_values <- denominator_data[i, -1]
        
        if (length(unique(num_values)) == 1 || length(unique(denom_values)) == 1 ||
            any(is.na(num_values)) || any(is.na(denom_values))) {
          p_values[i] <- NA
        } else {
          try({
            p_values[i] <- t.test(num_values, denom_values)$p.value
          }, silent = TRUE)
        }
      }
      
      data.frame(
        Compound_Name = numerator_data$Compound_Name,
        p_value = p_values,
        padj    = p.adjust(p_values, method = "BH"),
        row.names = NULL
      )
    })
    
    # UI for Carbon chain range
    output$carbon_range_ui <- renderUI({
      # Use all cleaned lipid names in this processed dataset
      nm_all <- compute_names_mapping(data[[1]])
      
      carbon_vals <- suppressWarnings(as.integer(nm_all$Carbon))
      carbon_vals <- carbon_vals[is.finite(carbon_vals)]
      if (length(carbon_vals) == 0) return(NULL)
      
      sliderInput(
        "carbon_range",
        "Carbon chain range:",
        min   = min(carbon_vals),
        max   = max(carbon_vals),
        value = c(min(carbon_vals), max(carbon_vals)),
        step  = 1
      )
    })
    
    # UI for Double-bond range
    output$double_range_ui <- renderUI({
      nm_all <- compute_names_mapping(data[[1]])
      
      double_vals <- suppressWarnings(as.integer(nm_all$Double))
      double_vals <- double_vals[is.finite(double_vals)]
      if (length(double_vals) == 0) return(NULL)
      
      sliderInput(
        "double_range",
        "Double bonds range:",
        min   = min(double_vals),
        max   = max(double_vals),
        value = c(min(double_vals), max(double_vals)),
        step  = 1
      )
    })
    
    
    
    
    
    # p-values for ALL lipids (no filtering); used only for drawing outlines
    reactiveP_value_all <- reactive({
      req(input$selected_group_for_numerator, input$selected_group_for_denominator)
      
      numerator_data   <- grouped_data_frames[[input$selected_group_for_numerator]]
      denominator_data <- grouped_data_frames[[input$selected_group_for_denominator]]
      
      if (nrow(numerator_data) == 0 || nrow(denominator_data) == 0) return(NULL)
      
      p_values <- numeric(nrow(numerator_data))
      
      for (i in 1:nrow(numerator_data)) {
        num_values   <- numerator_data[i, -1]
        denom_values <- denominator_data[i, -1]
        
        if (length(unique(num_values)) == 1 || length(unique(denom_values)) == 1 ||
            any(is.na(num_values)) || any(is.na(denom_values))) {
          p_values[i] <- NA
        } else {
          p_values[i] <- t.test(num_values, denom_values)$p.value
        }
      }
      
      data.frame(
        Compound_Name = numerator_data$Compound_Name,
        p_value = p_values,
        padj    = p.adjust(p_values, method = "BH"),
        row.names = NULL
      )
    })
    
    
    
    # Warring message if no data meets the threshold, but does not show if 'all' lipids are selected.
    output$filteredDataWarning <- renderUI({
      filtered_data <- reactiveFilteredData()
      
      if (selected_lipid_count() == 0) {
        # No lipids selected, do not show any message
        return(NULL)
      } else if (nrow(filtered_data) == 0 && selected_lipid_count() >= 1) {
        # No data meets the filtering criteria and lipids are selected
        div(
          style = "color: red; font-weight: bold;",
          "Warning: No data meets the filtering criteria."
        )
      } else {
        # Data exists or no lipids are selected, no warning needed
        NULL
      }
    })
    
    
    # Build overlay dataset with tile positions for significant lipids
    sig_points <- reactive({
      if (!isTRUE(input$outline_sig)) return(NULL)
      
      base_df <- reactiveFilteredData()       # the lipids currently shown in heatmap
      pv_all  <- reactiveP_value_all()
      if (is.null(base_df) || nrow(base_df) == 0 || is.null(pv_all)) return(NULL)
      
      # join to get both p_value and padj for the lipids actually plotted
      j <- merge(base_df["Compound_Name"], pv_all, by = "Compound_Name", all.x = TRUE)
      
      metric <- if (identical(input$sig_metric, "padj")) "padj" else "p_value"
      j <- j[!is.na(j[[metric]]) & j[[metric]] <= input$sig_threshold, , drop = FALSE]
      if (nrow(j) == 0) return(NULL)
      
      nm  <- map_lipid_names(x = j$Compound_Name)  # gives Class
      pos <- summarize_cd(j$Compound_Name)         # gives totalC, totalD
      
      df <- data.frame(
        Compound_Name = j$Compound_Name,
        Class  = nm$Class,
        totalC = pos[, "totalC"],
        totalD = pos[, "totalD"],
        stringsAsFactors = FALSE
      )
      
      df <- df[!is.na(df$totalC) & !is.na(df$totalD), , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
      df
    })
    
    
    
    
    ###### Creating the lipid heatmap plot
    heatmap_plot_data <- reactive({
      filtered_data <- reactiveFilteredData()
      req(nrow(filtered_data) > 0)
      
      num_group <- input$selected_group_for_numerator
      den_group <- input$selected_group_for_denominator
      
      ## --- MAIN TITLE ---------------------------------------------------------
      main_title <- if (!is.null(input$main_title) && nzchar(input$main_title)) {
        input$main_title
      } else {
        NULL
      }
      
      ## --- SUBTITLE: user text + custom group joiner -------------------------
      user_txt <- if (!is.null(input$sub_title)) trimws(input$sub_title) else ""
      join_txt <- if (!is.null(input$group_join)) trimws(input$group_join) else ""
      if (!nzchar(join_txt)) join_txt <- "vs"
      
      group_label <- if (!is.null(num_group) && nzchar(num_group) &&
                         !is.null(den_group) && nzchar(den_group)) {
        paste(num_group, join_txt, den_group)
      } else {
        ""
      }
      
      if (isTRUE(input$auto_subtitle_groups)) {
        if (nzchar(user_txt) && nzchar(group_label)) {
          subtitle_text <- paste(user_txt, group_label)
        } else if (!nzchar(user_txt) && nzchar(group_label)) {
          subtitle_text <- group_label
        } else if (nzchar(user_txt) && !nzchar(group_label)) {
          subtitle_text <- user_txt
        } else {
          subtitle_text <- NULL
        }
      } else {
        subtitle_text <- if (nzchar(user_txt)) user_txt else NULL
      }
      
      ## --- NAME MAPPING -------------------------------------------------------
      nm <- compute_names_mapping(filtered_data$Compound_Name)
      
      plot_df <- dplyr::left_join(
        filtered_data,
        nm[, c("Compound_Name", "Class", "Carbon", "Double")],
        by = "Compound_Name"
      )
      plot_df$Class  <- nm$Class[match(plot_df$Compound_Name, nm$Compound_Name)]
      plot_df$Carbon <- as.integer(plot_df$Carbon)
      plot_df$Double <- as.integer(plot_df$Double)
      
      # Facet variable (handles isoforms: TG_1(58:8) -> TG_1)
      plot_df$FacetClass <- sub("\\s*\\(.*$", "", plot_df$Compound_Name)
      plot_df$FacetClass <- droplevels(factor(plot_df$FacetClass))
      
      ## --- facet-wise padding to avoid tight axes -----------------------------
      min_span_x <- 3
      min_span_y <- 3
      
      pad_tbl <- plot_df %>%
        dplyr::group_by(FacetClass) %>%
        dplyr::summarise(
          xmin = min(Carbon, na.rm = TRUE),
          xmax = max(Carbon, na.rm = TRUE),
          ymin = min(Double, na.rm = TRUE),
          ymax = max(Double, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          xmin = floor(xmin - 0.5), xmax = ceiling(xmax + 0.5),
          ymin = floor(ymin - 0.5), ymax = ceiling(ymax + 0.5),
          span_x = pmax(xmax - xmin, min_span_x),
          span_y = pmax(ymax - ymin, min_span_y),
          cx = (xmin + xmax) / 2, cy = (ymin + ymax) / 2,
          xmin = floor(cx - span_x / 2), xmax = ceiling(cx + span_x / 2),
          ymin = floor(cy - span_y / 2), ymax = ceiling(cy + span_y / 2)
        )
      
      pad_df <- rbind(
        data.frame(FacetClass = pad_tbl$FacetClass,
                   Carbon = pad_tbl$xmin, Double = pad_tbl$ymin),
        data.frame(FacetClass = pad_tbl$FacetClass,
                   Carbon = pad_tbl$xmax, Double = pad_tbl$ymax)
      )
      
      ## --- facet layout -------------------------------------------------------
      ## --- DYNAMIC CALCULATION ------------------------------------------------
      # 1. Count how many unique facets (lipid classes) we have
      n_facets <- length(unique(plot_df$FacetClass))
      
      # 2. Determine Columns
      # If Manual: use input. If Auto: default to 3.
      if (input$layout_mode == "custom" && !is.null(input$facet_cols)) {
        ncol_facets <- as.integer(input$facet_cols)
      } else {
        ncol_facets <- 3L
      }
      
      # 3. Calculate necessary Rows
      nrow_facets <- ceiling(n_facets / ncol_facets)
      
      ## --- WIDTH/HEIGHT LOGIC -------------------------------------------------
      if (input$layout_mode == "auto") {
        # --- AUTO MODE ---
        # Calculate dynamic height: 
        # Assign ~300px per row of charts + 150px buffer for titles/legends
        # Assign ~350px per column + buffer
        
        calc_height <- (nrow_facets * 300) + 150
        calc_width  <- (ncol_facets * 350) + 100
        
        # Set limits so it doesn't get too small or insanely large automatically
        plot_height_px <- max(600, min(12000, calc_height))
        plot_width_px  <- max(1000, min(5000, calc_width))
        
      } else {
        # --- CUSTOM MODE ---
        # Use user inputs, with fallbacks
        plot_height_px <- if (!is.null(input$plot_height_px) && 
                              !is.na(input$plot_height_px) && 
                              input$plot_height_px > 0) {
          input$plot_height_px
        } else {
          1000
        }
        
        plot_width_px <- if (!is.null(input$plot_width_px) && 
                             !is.null(input$plot_width_px) && 
                             input$plot_width_px > 0) {
          input$plot_width_px
        } else {
          1100
        }
      }
      
      # clamp (keep existing safety limits)
      min_h <- 400; max_h <- 15000 # Increased max_h slightly for big datasets
      min_w <- 400; max_w <- 6000
      plot_height_px <- max(min_h, min(max_h, plot_height_px))
      plot_width_px <- max(min_w, min(max_w, plot_width_px))
      
      ## --- Color limits / clamping --------------------------------------------
      logFC_range <- range(plot_df$logFC, na.rm = TRUE)
      if (input$selected_logfc_sclae_bar == "Manual") {
        lim <- input$logFC_scale_manual
        fill.limits <- c(-lim, lim)
      } else {
        fill.limits <- c(min(logFC_range), max(logFC_range))
      }
      
      plot_df$fill_clamped <- pmax(pmin(plot_df$logFC, fill.limits[2]), fill.limits[1])
      plot_df$extreme_label <- ifelse(
        plot_df$logFC > fill.limits[2], "+",
        ifelse(plot_df$logFC < fill.limits[1], "-", "")
      )
      
      ## --- tick helper --------------------------------------------------------
      int_breaks <- function(n = 6) {
        function(lims) {
          rng <- lims[2] - lims[1]
          if (!is.finite(rng)) return(lims)
          step <- max(1, ceiling(rng / n))
          seq(floor(lims[1]), ceiling(lims[2]), by = step)
        }
      }
      
      ## --- significance overlay -----------------------------------------------
      sig_df <- NULL
      if (isTRUE(input$outline_sig)) {
        metric_col <- if (identical(input$sig_metric, "padj")) "padj" else "p_value"
        cutoff <- input$sig_threshold
        if (is.null(cutoff) || is.na(cutoff)) cutoff <- 0.05
        cutoff <- max(0, min(1, as.numeric(cutoff)))
        
        if (metric_col %in% names(plot_df)) {
          sig_df <- dplyr::filter(
            plot_df,
            !is.na(.data[[metric_col]]) & .data[[metric_col]] <= cutoff
          )
          if (nrow(sig_df) == 0) sig_df <- NULL
        }
        
        if (!is.null(sig_df) && identical(input$sig_shape, "circle")) {
          mval <- pmax(sig_df[[metric_col]], .Machine$double.xmin)
          sig_df$size_raw <- -log10(mval)
          min_mm <- 2.5; max_mm <- 9
          rng <- range(sig_df$size_raw[is.finite(sig_df$size_raw)], na.rm = TRUE)
          sig_df$size_bubble <- if (is.finite(rng[1]) && is.finite(rng[2]) && diff(rng) > 0) {
            min_mm + (sig_df$size_raw - rng[1]) * (max_mm - min_mm) / diff(rng)
          } else {
            (min_mm + max_mm) / 2
          }
        }
      }
      
      ## --- tile border settings (already installed) ---------------------------
      border_color <- if (isTRUE(input$tile_border_show) &&
                          !is.null(input$tile_border_color) &&
                          !is.na(input$tile_border_size) &&
                          input$tile_border_size > 0) {
        input$tile_border_color
      } else {
        NA
      }
      
      border_size <- if (isTRUE(input$tile_border_show) &&
                         !is.na(input$tile_border_size) &&
                         input$tile_border_size > 0) {
        input$tile_border_size
      } else {
        0
      }
      
      ## --- LEGEND SIDE + TITLE ALIGNMENT (NEW) --------------------------------
      legend_side <- if (!is.null(input$legend_side)) tolower(input$legend_side) else "top"
      if (legend_side %in% c("none", "hide")) {
        legend_pos_value <- "none"
      } else if (legend_side %in% c("top", "bottom", "left", "right")) {
        legend_pos_value <- legend_side
      } else {
        legend_pos_value <- "top"
      }
      
      colorbar_direction <- if (legend_pos_value %in% c("left", "right")) "vertical" else "horizontal"
      label_position     <- if (colorbar_direction == "horizontal") "bottom" else "right"
      
      title_hjust    <- if (!is.null(input$title_hjust))    input$title_hjust    else 0.5
      subtitle_hjust <- if (!is.null(input$subtitle_hjust)) input$subtitle_hjust else 0.5
      
      # NEW: base font family for all text
      base_family <- if (!is.null(input$font_family) &&
                         nzchar(input$font_family) &&
                         input$font_family != "Default") {
        input$font_family
      } else {
        NULL  # use ggplot2 / R default
      }
      
      
      # --- Compute bar dimensions depending on orientation ---------------------
      bw <- input$barwidth
      bh <- input$barheight
      
      if (colorbar_direction == "vertical") {
        # When the bar is vertical (legend on left/right), flip width/height
        tmp <- bw
        bw  <- bh
        bh  <- tmp
      }
      
      ## --- base plot -----------------------------------------------------------
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(Carbon, Double)) +
        ggplot2::geom_blank(data = pad_df,
                            ggplot2::aes(x = Carbon, y = Double)) +
        ggplot2::geom_tile(
          ggplot2::aes(fill = fill_clamped),
          width  = 0.95,
          height = 0.95,
          color  = border_color,
          linewidth = border_size
        ) +
        ggplot2::geom_text(ggplot2::aes(label = extreme_label), na.rm = TRUE) +
        ggplot2::scale_x_continuous(
          breaks = int_breaks(6),
          expand = ggplot2::expansion(mult = c(0.01, 0.01))
        ) +
        ggplot2::scale_y_continuous(
          breaks = int_breaks(6),
          expand = ggplot2::expansion(mult = c(0.01, 0.01))
        ) +
        ggplot2::scale_fill_gradient2(
          low    = input$low_color,
          mid    = input$mid_color,
          high   = input$high_color,
          limits = fill.limits,
          space  = "Lab",
          name   = expression(log[2](FC)),
          guide  = if (identical(legend_pos_value, "none")) "none" else
            ggplot2::guide_colorbar(
              direction      = colorbar_direction,  # horizontal vs vertical
              title.position = "top",
              title.hjust    = 0.5,
              label.position = label_position,
              barwidth       = bw,
              barheight      = bh
            )
        ) +
        ggplot2::facet_wrap(
          ~ FacetClass,
          scales = "free",
          ncol   = ncol_facets,
          drop   = TRUE
        ) +
        ggplot2::theme(
          # Base font for all text
          text = ggplot2::element_text(family = base_family),
          
          panel.background = ggplot2::element_rect(
            fill = input$panel_bg_color, color = "white"
          ),
          strip.background = ggplot2::element_rect(
            fill = input$strip_bg_color, color = "white"
          ),
          strip.text = ggplot2::element_text(
            color = input$strip_text_color,
            face  = "bold",
            size  = input$strip_text_size
          ),
          axis.text.x = ggplot2::element_text(
            angle = input$axis_text_x_angle,
            size  = input$axis_text_x_size
          ),
          axis.text.y  = ggplot2::element_text(size = input$axis_text_y_size),
          axis.title   = ggplot2::element_text(size = input$axis_title_size),
          legend.title = ggplot2::element_text(size = input$legend_title_size),
          legend.text  = ggplot2::element_text(size = input$legend_text_size),
          axis.title.x = ggplot2::element_text(size = input$axis_title_x_size),
          axis.title.y = ggplot2::element_text(size = input$axis_title_y_size),
          plot.title   = ggplot2::element_text(
            size  = input$plot_title_size,
            face  = "bold",
            hjust = title_hjust
          ),
          plot.subtitle = ggplot2::element_text(
            size  = input$plot_subtitle_size,
            face  = if (isTRUE(input$plot_subtitle_bold)) "bold" else "plain",
            hjust = subtitle_hjust
          ),
          
          # --- NEW: spacing & margins in mm ---
          panel.spacing = grid::unit(
            if (!is.null(input$panel_spacing_mm)) input$panel_spacing_mm else 2,
            "mm"
          ),
          plot.margin = ggplot2::margin(
            t = if (!is.null(input$margin_top_mm))    input$margin_top_mm    else 5,
            r = if (!is.null(input$margin_right_mm))  input$margin_right_mm  else 5,
            b = if (!is.null(input$margin_bottom_mm)) input$margin_bottom_mm else 5,
            l = if (!is.null(input$margin_left_mm))   input$margin_left_mm   else 5,
            unit = "mm"
          ),
          
          
          legend.position = legend_pos_value
        ) +
        
        ggplot2::labs(
          x        = input$x_axis_label,
          y        = input$y_axis_label,
          title    = main_title,
          subtitle = subtitle_text
        )
      
      ## --- significance overlays ----------------------------------------------
      if (!is.null(sig_df)) {
        if (identical(input$sig_shape, "circle")) {
          p <- p +
            ggplot2::geom_point(
              data = sig_df,
              ggplot2::aes(x = Carbon, y = Double, size = size_bubble),
              shape  = 21,
              fill   = NA,
              colour = if (!is.null(input$sig_outline_color)) input$sig_outline_color else "navy",
              stroke = if (!is.null(input$sig_outline_size)) input$sig_outline_size else 0.7,
              show.legend = FALSE
            ) +
            ggplot2::scale_size_identity()
        } else {
          p <- p +
            ggplot2::geom_tile(
              data = sig_df,
              ggplot2::aes(x = Carbon, y = Double),
              fill      = NA,
              color     = if (!is.null(input$sig_outline_color)) input$sig_outline_color else "navy",
              linewidth = if (!is.null(input$sig_outline_size)) input$sig_outline_size else 0.7,
              width     = 0.95,
              height    = 0.95,
              show.legend = FALSE
            )
        }
      }
      
      ## --- grid on/off --------------------------------------------------------
      if (!isTRUE(input$show_grid)) {
        # No grid lines at all
        p <- p + ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()
        )
      } else {
        # Show major grid lines in user-selected color, no minor grid
        p <- p + ggplot2::theme(
          panel.grid.major = ggplot2::element_line(
            color = input$grid_color,
            linewidth = 0.3
          ),
          panel.grid.minor = ggplot2::element_blank()
        )
      }
      
      
      list(plot = p, height = plot_height_px, width = plot_width_px)
    })
    
    
    
    
    
    
    # Render the heatmap plot with dynamic width and height (in px)
    output$heatmapPlot <- renderPlot({
      heatmap_plot_data()$plot
    },
    height = function() {
      heatmap_plot_data()$height
    },
    width = function() {
      heatmap_plot_data()$width
    })
    
    
    
    
    
    # Layout for visualization (uses the same dynamic height; width on the plot itself)
    output$visualization_ui <- renderUI({
      plot_info <- heatmap_plot_data()
      plot_h <- plot_info$height  # px
      plot_w <- plot_info$width   # px
      
      if (input$split_screen) {
        # Heatmap + table side by side, same vertical space
        fluidRow(
          column(
            width = 6,
            div(
              style = paste0("width: 100%; height: ", plot_h, "px; overflow-y: auto;"),
              withSpinner(
                plotOutput(
                  "heatmapPlot",
                  width  = paste0(plot_w, "px"),
                  height = paste0(plot_h, "px")
                )
              )
            )
          ),
          column(
            width = 6,
            div(
              style = paste0("width: 100%; height: ", plot_h, "px; overflow-y: auto;"),
              withSpinner(dataTableOutput("pValueTable_2"))
            )
          )
        )
      } else {
        # Only heatmap
        div(
          style = paste0("width: 100%; height: ", plot_h, "px; overflow-y: auto;"),
          withSpinner(
            plotOutput(
              "heatmapPlot",
              width  = paste0(plot_w, "px"),
              height = paste0(plot_h, "px"),
              click = "heatmap_click"
            )
          )
        )
      }
    })
    
    
    
    
    # ---- Color palette presets (CB-friendly) -------------------------------
    observeEvent(input$palette_preset, {
      pal <- input$palette_preset
      if (is.null(pal) || pal == "manual") return()
      
      # Each branch overwrites low / mid / high logFC colors
      if (pal == "default") {
        # Your current default
        updateColourInput(session, "low_color",  value = "#4575b4")
        updateColourInput(session, "mid_color",  value = "#FFFFFF")
        updateColourInput(session, "high_color", value = "#d73027")
        
      } else if (pal == "blue_orange") {
        # Okabe–Ito blue vs vermillion (color-blind friendly)
        updateColourInput(session, "low_color",  value = "#0072B2")  # blue
        updateColourInput(session, "mid_color",  value = "#FFFFFF")
        updateColourInput(session, "high_color", value = "#D55E00")  # vermillion
        
      } else if (pal == "purple_green") {
        # Okabe–Ito purple vs green (color-blind friendly)
        updateColourInput(session, "low_color",  value = "#CC79A7")  # purple
        updateColourInput(session, "mid_color",  value = "#FFFFFF")
        updateColourInput(session, "high_color", value = "#009E73")  # green
        
      } else if (pal == "viridis_like") {
        # Viridis-inspired sequential (good for “all positive” heatmaps too)
        updateColourInput(session, "low_color",  value = "#440154")  # dark purple
        updateColourInput(session, "mid_color",  value = "#21918C")  # teal
        updateColourInput(session, "high_color", value = "#FDE725")  # yellow
      }
    })
    
    
    
    #### ---- Reset buttons for Plot settings tabs ---- ####
    
    # Colors tab
    observeEvent(input$reset_colors, {
      updateSelectInput(session, "palette_preset", selected = "default")
      
      updateColourInput(session, "low_color",  value = "#4575b4")
      updateColourInput(session, "mid_color",  value = "#FFFFFF")
      updateColourInput(session, "high_color", value = "#d73027")
      
      updateColourInput(session, "panel_bg_color", value = "#D3D3D3")
      updateColourInput(session, "strip_bg_color", value = "#3483d1")
      updateColourInput(session, "strip_text_color", value = "black")
      
      updateCheckboxInput(session, "tile_border_show", value = TRUE)
      updateColourInput(session, "tile_border_color", value = "white")
      updateNumericInput(session, "tile_border_size", value = 0.2)
    })
    
    
    # Text Sizes tab
    observeEvent(input$reset_text_sizes, {
      updateNumericInput(session, "strip_text_size",    value = 16)
      updateNumericInput(session, "axis_text_x_size",   value = 12)
      updateNumericInput(session, "axis_text_y_size",   value = 12)
      updateNumericInput(session, "axis_title_size",    value = 20)
      updateNumericInput(session, "legend_title_size",  value = 18)
      updateNumericInput(session, "legend_text_size",   value = 14)
      updateNumericInput(session, "axis_title_x_size",  value = 20)
      updateNumericInput(session, "axis_title_y_size",  value = 20)
      updateNumericInput(session, "plot_title_size",    value = 22)
      updateNumericInput(session, "plot_subtitle_size", value = 18)
      updateCheckboxInput(session, "plot_subtitle_bold", value = FALSE)
    })
    
    # Axis Settings tab
    observeEvent(input$reset_axis, {
      updateNumericInput(session, "axis_text_x_angle", value = 90)
    })
    
    # Legend Settings tab
    observeEvent(input$reset_legend, {
      updateNumericInput(session, "barwidth",  value = 15)
      updateNumericInput(session, "barheight", value = 1)
    })
    
    # Labels tab
    observeEvent(input$reset_labels, {
      updateTextInput(session, "x_axis_label",
                      value = "Number of fatty-acid carbon atoms")
      updateTextInput(session, "y_axis_label",
                      value = "Number of fatty-acid double bonds")
      updateTextInput(session, "main_title",
                      value = "Lipid Fold Changes by Class")
      updateTextInput(session, "sub_title",
                      value = "Groups being compared is")
      updateTextInput(session, "group_join",
                      value = "vs")
      updateCheckboxInput(session, "auto_subtitle_groups", value = TRUE)
    })
    
    # Layout tab
    observeEvent(input$reset_layout, {
      updateNumericInput(session, "facet_cols",      value = 3)
      updateRadioButtons(session, "layout_mode", selected = "auto")
      updateNumericInput(session, "plot_width_px",   value = 1100)
      updateNumericInput(session, "plot_height_px",  value = 3000)
      updateNumericInput(session, "panel_spacing_mm",  value = 2)
      updateNumericInput(session, "margin_top_mm",     value = 5)
      updateNumericInput(session, "margin_right_mm",   value = 5)
      updateNumericInput(session, "margin_bottom_mm",  value = 5)
      updateNumericInput(session, "margin_left_mm",    value = 5)
    })
    
    
    
    
    observeEvent(input$heatmap_click, {
      # 1. Get the clicked coordinates
      click <- input$heatmap_click
      
      # Return early if clicked outside or no plot exists
      if (is.null(click$x) || is.null(click$y)) return()
      
      # 2. Fetch the fully processed data DIRECTLY from the plot object
      # This ensures we have 'Carbon', 'Double', and 'FacetClass' columns
      plot_result <- heatmap_plot_data()
      req(plot_result$plot$data)
      plot_df <- plot_result$plot$data
      
      # 3. Match the click to the data
      # Round coordinates because Carbon/Double are integers on the axis
      clicked_carbon <- round(click$x)
      clicked_double <- round(click$y)
      
      # Handle Facets: ggplot sends the facet value in 'panelvar1'
      clicked_facet  <- click$panelvar1
      
      # Find the specific lipid
      selected_lipid <- plot_df %>%
        dplyr::filter(
          Carbon == clicked_carbon,
          Double == clicked_double,
          # Only filter by facet if the click contains facet info (it usually does for faceted plots)
          if (!is.null(clicked_facet)) FacetClass == clicked_facet else TRUE
        )
      
      # If a unique lipid is found, show the modal
      if (nrow(selected_lipid) == 1) {
        lipid_name <- selected_lipid$Compound_Name
        
        # 4. Retrieve Raw Data for this lipid from the grouped data frames
        num_group <- input$selected_group_for_numerator
        den_group <- input$selected_group_for_denominator
        
        # Helper to extract values for a specific lipid from a group dataframe
        get_values <- function(group_name, lipid) {
          # Access the reactive value containing the grouped data
          df <- values$grouped_data_frames[[group_name]]
          
          # Find row with this lipid
          row_idx <- which(df$Compound_Name == lipid)
          
          if (length(row_idx) > 0) {
            # Extract numeric values (excluding Compound_Name at col 1)
            val <- as.numeric(df[row_idx, -1]) 
            return(data.frame(Group = group_name, Abundance = val))
          }
          return(NULL)
        }
        
        df_num <- get_values(num_group, lipid_name)
        df_den <- get_values(den_group, lipid_name)
        plot_data <- rbind(df_num, df_den)
        
        # 5. Render the Modal
        showModal(modalDialog(
          title = paste("Detail View:", lipid_name),
          renderPlot({
            req(plot_data)
            ggplot(plot_data, aes(x = Group, y = Abundance, fill = Group)) +
              geom_boxplot(alpha = 0.6, outlier.shape = NA) +
              geom_jitter(width = 0.2, size = 2) +
              theme_minimal(base_size = 14) +
              labs(
                title = paste0("LogFC: ", round(selected_lipid$logFC, 3), 
                               " | p-adj: ", formatC(selected_lipid$padj, format = "e", digits = 2)),
                subtitle = paste("Carbon:", selected_lipid$Carbon, "| Double Bonds:", selected_lipid$Double),
                y = "Abundance (Intensity)"
              ) +
              scale_fill_manual(values = c(input$high_color, input$low_color)) 
          }),
          size = "m",
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }
    })
    
    
    
    
    # Reactive expression to calculate lipid group summary and total lipid count
    lipid_summary <- reactive({
      lipid_group_df <- as.data.frame(table(group_lipids_by_group(data)$group)) 
      colnames(lipid_group_df) <- c("Lipid group", "Count")
      
      # Calculate and add percentage of total lipids
      lipid_group_df$`Percentage of Total` <- round((lipid_group_df$Count / sum(lipid_group_df$Count)) * 100, 2)
      
      return(lipid_group_df)  # Return the dataframe for use in other parts
    })
    
    # Function to generate lipid group summary and percentage for the table
    output$lipid_group_count <- DT::renderDataTable({
      lipid_group_df <- lipid_summary()  # Use the reactive expression
      
      # Render the table
      DT::datatable(lipid_group_df, options = list(pageLength = 5, autoWidth = TRUE))
    })
    
    # Render total lipids as a text output
    output$lipid_total <- renderText({
      lipid_group_df <- lipid_summary()  # Use the reactive expression
      lipidsum <- sum(lipid_group_df$Count)
      paste("Total lipid count in table (After data cleaning, and are not affect by threshold set in 'Lipid Visualization'):", lipidsum)
    })
    
    
    # Assume merged_data_info has columns: Compound_Name, merged_molecules, count
    # returned by merged_info_function(data)
    
    # render table with "Compound_Name", "Original.annotation", "logFC", "p_value", "padj"
    output$pValueTable <- renderDataTable({
      filtered_data <- reactiveFilteredData()
      req(data)  # We must have data
      
      # Ensure 'Compound_Name' exists in 'data'
      if (!"Compound_Name" %in% colnames(data)) {
        # FIX: Use [[1]] to get a vector, avoiding tibble/dataframe nested column issues
        data$Compound_Name <- data[[1]]  
      }
      
      # Ensure 'Compound_Name' exists in 'filtered_data'
      if (!"Compound_Name" %in% colnames(filtered_data)) {
        stop("'Compound_Name' column not found in 'filtered_data'")
      }
      
      # Check if 'Original annotation' exists in 'data' and rename
      if ("Original annotation" %in% colnames(data)) {
        colnames(data)[colnames(data) == "Original annotation"] <- "Original.annotation"
      }
      
      # Merge filtered_data with Original.annotation if it exists
      if ("Original.annotation" %in% colnames(data)) {
        merged_data <- merge(
          filtered_data,
          data[, c("Compound_Name", "Original.annotation")],
          by = "Compound_Name",
          all.x = TRUE
        )
      } else {
        merged_data <- filtered_data
      }
      
      # Merge in the merged_data_info if available
      if (!is.null(merged_data_info) && "Compound_Name" %in% colnames(merged_data_info)) {
        merged_data <- merge(
          merged_data,
          merged_data_info,
          by = "Compound_Name",
          all.x = TRUE
        )
      }
      
      # Determine which columns to show:
      columns_to_show <- c("Compound_Name", "logFC", "p_value", "padj")
      
      if ("Original.annotation" %in% colnames(merged_data)) {
        columns_to_show <- c("Compound_Name", "Original.annotation", columns_to_show)
      }
      
      if ("merged_molecules" %in% colnames(merged_data)) {
        # Include merged_molecules and count columns
        columns_to_show <- c("Compound_Name", "Original.annotation", "merged_molecules", "count", "logFC", "p_value", "padj")
      }
      
      # Create the final table to show
      dataTableToShow <- merged_data[, intersect(columns_to_show, colnames(merged_data)), drop = FALSE]
      
      # Round numeric columns if they exist
      numeric_cols <- intersect(c("logFC", "p_value", "padj"), colnames(dataTableToShow))
      dataTableToShow[numeric_cols] <- lapply(dataTableToShow[numeric_cols], function(x) round(x, 5))
      
      # Render the DataTable
      datatable(dataTableToShow, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    
    
    # This table is smaler and does not contain the 'Original annotation' column, is used to split the screen with table and heatmap
    output$pValueTable_2 <- renderDataTable({
      filtered_data <- reactiveFilteredData()
      
      
      dataTableToShow <- filtered_data[, c("Compound_Name", "logFC", "p_value", "padj")]
      
      # Round 'logFC' and 'p_value' to the desired number of decimal places
      dataTableToShow$logFC <- round(dataTableToShow$logFC, 5)      # 5 decimal places for logFC
      dataTableToShow$p_value <- round(dataTableToShow$p_value, 5)  # 5 decimal places for p-value
      dataTableToShow$padj <- round(dataTableToShow$padj, 5)  # 5 decimal places for p-value
      
      
      # Render the selected data in a DataTable
      datatable(dataTableToShow, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    
    base_dpi <- 72  # dpi used by Shiny's renderPlot (approx.)
    
    plot_dimensions <- reactive({
      plot_info <- heatmap_plot_data()
      req(plot_info$width, plot_info$height)  # these are in *pixels*
      
      # Physical size in inches as drawn in the app
      width_in  <- plot_info$width  / base_dpi
      height_in <- plot_info$height / base_dpi
      
      # Optional: clamp to something reasonable
      width_in  <- max(3, min(40, width_in))
      height_in <- max(3, min(40, height_in))
      
      list(width_in = width_in, height_in = height_in)
    })
    
    
    
    
    # Observe the action button to open the modal dialog
    observeEvent(input$download_heatmap_btn, {
      showModal(modalDialog(
        title = "Download Heatmap Image",
        
        # File Name Input
        textInput("modal_file_name", "Enter File Name", value = paste0("lipid_heatmap_", Sys.Date())),
        
        # Format Selection
        selectInput("modal_image_format", "Select Image Format", choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
        
        # --- UI CHANGE START ---
        # 1. Add 'max = 600' so the spinner stops there.
        numericInput("modal_image_dpi", "Image Resolution (DPI)", value = 300, min = 72, max = 600, step = 72),
        
        # 2. Add a helpful note so the user sees the difference immediately.
        tags$p(style = "color: #888; font-size: 12px; margin-top: -10px;", 
               "Note: Maximum resolution for PNG/JPEG is 600 DPI."),
        # --- UI CHANGE END ---
        
        footer = tagList(
          modalButton("Cancel"),
          downloadButton("modal_download_heatmap", "Download")
        )
      ))
    })
    
    # ---- FIX: STRICT INPUT ENFORCEMENT ----
    # This observer watches what the user types. 
    # If they go above 600 for PNG/JPEG, it forces the number back down.
    
    observeEvent(c(input$modal_image_dpi, input$modal_image_format), {
      # Only run if the modal is open and input exists
      req(input$modal_image_dpi)
      
      # Define the safety limit
      max_safe_dpi <- 600
      
      # Check if we need to enforce the limit (Raster images only)
      is_vector <- tolower(input$modal_image_format) %in% c("pdf", "svg")
      
      if (!is_vector && input$modal_image_dpi > max_safe_dpi) {
        
        # 1. Snap the value back to 600 immediately
        updateNumericInput(session, "modal_image_dpi", value = max_safe_dpi)
        
        # 2. Show a polite, transient warning so they know why it changed
        showNotification(
          paste("Resolution capped at", max_safe_dpi, "DPI for this format."), 
          type = "warning", 
          duration = 3
        )
      }
    })
    # ---------------------------------------
    
    
    
    
    
    
    
    output$modal_download_heatmap <- downloadHandler(
      filename = function() {
        # Retrieve the name from the modal input
        base_name <- if (nzchar(input$modal_file_name)) input$modal_file_name else "heatmap_plot"
        
        # Combine user name with the selected extension
        paste0(base_name, ".", input$modal_image_format)
      },
      content = function(file) {
        hp   <- heatmap_plot_data()
        plot <- hp$plot
        
        dims <- plot_dimensions()
        
        
        # 1. Get the user's desired DPI
        user_dpi <- if (!is.null(input$modal_image_dpi) &&
                        !is.na(input$modal_image_dpi) &&
                        input$modal_image_dpi > 0) {
          input$modal_image_dpi
        } else {
          300
        }
        
        # 2. Apply limits based on file type
        # If it is a Raster image (PNG, JPEG), we cap it at 600 DPI to prevent server crashes.
        # If it is a Vector image (PDF, SVG), we allow the user's input (or default to 300) 
        # because vectors are math-based and won't consume RAM like pixels do.
        
        if (tolower(input$modal_image_format) %in% c("pdf", "svg")) {
          dpi <- user_dpi 
        } else {
          dpi <- min(user_dpi, 600) # <--- THE SAFETY CAP
        }
        
        ggsave(
          filename = file,
          plot     = plot,
          device   = input$modal_image_format,
          width    = dims$width_in,
          height   = dims$height_in,
          units    = "in",
          dpi      = dpi
        )
        
        removeModal() # Closes modal after download starts
      }
    )      
    
    
    
    
    # ---- Supplementary table (logFC + p + padj + annotation) ---------------
    supplementary_table <- reactive({
      filtered_data <- reactiveFilteredData()
      req(data)  # original/merged data available
      
      # Ensure 'Compound_Name' exists in 'data'
      if (!"Compound_Name" %in% colnames(data)) {
        data$Compound_Name <- data[, 1]  # assume first column is name
      }
      
      # Make sure original annotation column name is consistent
      if ("Original annotation" %in% colnames(data)) {
        colnames(data)[colnames(data) == "Original annotation"] <- "Original.annotation"
      }
      
      # Start from filtered data
      merged_data <- filtered_data
      
      # Add Original.annotation if present
      if ("Original.annotation" %in% colnames(data)) {
        merged_data <- merge(
          merged_data,
          data[, c("Compound_Name", "Original.annotation")],
          by = "Compound_Name",
          all.x = TRUE
        )
      }
      
      # Add merged_data_info (for merged isoforms) if available
      if (!is.null(merged_data_info) &&
          "Compound_Name" %in% colnames(merged_data_info)) {
        merged_data <- merge(
          merged_data,
          merged_data_info,
          by = "Compound_Name",
          all.x = TRUE
        )
      }
      
      # Decide which columns to show
      columns_to_show <- c("Compound_Name", "logFC", "p_value", "padj")
      
      if ("Original.annotation" %in% colnames(merged_data)) {
        columns_to_show <- c("Compound_Name", "Original.annotation", columns_to_show)
      }
      
      if ("merged_molecules" %in% colnames(merged_data)) {
        columns_to_show <- c("Compound_Name",
                             "Original.annotation",
                             "merged_molecules",
                             "count",
                             "logFC", "p_value", "padj")
      }
      
      # Keep only available columns
      dataTableToShow <- merged_data[, intersect(columns_to_show, colnames(merged_data)), drop = FALSE]
      
      # Round numeric columns
      numeric_cols <- intersect(c("logFC", "p_value", "padj"), colnames(dataTableToShow))
      dataTableToShow[numeric_cols] <- lapply(
        dataTableToShow[numeric_cols],
        function(x) round(x, 5)
      )
      
      dataTableToShow
    })
    
    # On-screen table
    output$pValueTable <- DT::renderDataTable({
      DT::datatable(
        supplementary_table(),
        options = list(pageLength = 10, scrollX = TRUE)
      )
    })
    
    
    
    # ---- Download supplementary table as CSV -------------------------------
    output$download_pvalues_csv <- downloadHandler(
      filename = function() {
        paste0("LipidHeatmap_supplementary_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(supplementary_table(), file, row.names = FALSE)
      }
    )
    
    # ---- Download supplementary table as Excel (.xlsx) ---------------------
    output$download_pvalues_xlsx <- downloadHandler(
      filename = function() {
        paste0("LipidHeatmap_supplementary_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        # Needs: library(writexl)
        writexl::write_xlsx(supplementary_table(), path = file)
      }
    )
    
    
    
    
    
    
    #### upodate this with new features 
    # ---- Bundle heatmap style/layout settings into a list -------------------
    heatmap_settings <- reactive({
      list(
        # Thresholds & selection
        selected_lipid          = input$selected_lipid,
        logFC_input             = input$logFC_input,
        p_value_max             = input$p_value_max,
        p_value_adj             = input$p_value_adj,
        min_lipids_per_class    = input$min_lipids_per_class,
        
        # logFC scale & display
        selected_logfc_sclae_bar = input$selected_logfc_sclae_bar,
        logFC_scale_manual       = input$logFC_scale_manual,
        split_screen             = input$split_screen,
        show_grid                = input$show_grid,
        
        # significance overlay
        outline_sig        = input$outline_sig,
        sig_shape          = input$sig_shape,
        sig_metric         = input$sig_metric,
        sig_threshold      = input$sig_threshold,
        sig_outline_color  = input$sig_outline_color,
        sig_outline_size   = input$sig_outline_size,
        
        # color palette
        palette_preset     = input$palette_preset,
        low_color          = input$low_color,
        mid_color          = input$mid_color,
        high_color         = input$high_color,
        panel_bg_color     = input$panel_bg_color,
        strip_bg_color     = input$strip_bg_color,
        strip_text_color   = input$strip_text_color,
        tile_border_show   = input$tile_border_show,
        tile_border_color  = input$tile_border_color,
        tile_border_size   = input$tile_border_size,
        
        # text sizes & fonts
        strip_text_size    = input$strip_text_size,
        axis_text_x_size   = input$axis_text_x_size,
        axis_text_y_size   = input$axis_text_y_size,
        axis_title_size    = input$axis_title_size,
        legend_title_size  = input$legend_title_size,
        legend_text_size   = input$legend_text_size,
        axis_title_x_size  = input$axis_title_x_size,
        axis_title_y_size  = input$axis_title_y_size,
        plot_title_size    = input$plot_title_size,
        plot_subtitle_size = input$plot_subtitle_size,
        plot_subtitle_bold = input$plot_subtitle_bold,
        font_family        = input$font_family,
        
        # axis / legend
        axis_text_x_angle  = input$axis_text_x_angle,
        barwidth           = input$barwidth,
        barheight          = input$barheight,
        legend_side        = input$legend_side,
        
        # labels & titles
        x_axis_label       = input$x_axis_label,
        y_axis_label       = input$y_axis_label,
        main_title         = input$main_title,
        sub_title          = input$sub_title,
        group_join         = input$group_join,
        auto_subtitle_groups = input$auto_subtitle_groups,
        title_hjust        = input$title_hjust,
        subtitle_hjust     = input$subtitle_hjust,
        
        # layout
        facet_cols         = input$facet_cols,
        plot_width_px      = input$plot_width_px,
        plot_height_px     = input$plot_height_px,
        panel_spacing_mm   = input$panel_spacing_mm,
        margin_top_mm      = input$margin_top_mm,
        margin_right_mm    = input$margin_right_mm,
        margin_bottom_mm   = input$margin_bottom_mm,
        margin_left_mm     = input$margin_left_mm
      )
    })
    
    
    # ---- Export settings as .rds -------------------------------------------
    output$download_heatmap_settings <- downloadHandler(
      filename = function() {
        paste0("LipidHeatmap_settings_", Sys.Date(), ".rds")
      },
      content = function(file) {
        saveRDS(heatmap_settings(), file = file)
      }
    )
    
    
    
    # ---- Import settings from .rds and apply to inputs ---------------------
    observeEvent(input$upload_heatmap_settings, {
      req(input$upload_heatmap_settings$datapath)
      s <- readRDS(input$upload_heatmap_settings$datapath)
      if (!is.list(s)) return()
      
      # Helper: safely update if value exists in file
      safe_val <- function(name) {
        if (!is.null(s[[name]])) s[[name]] else NULL
      }
      
      # Thresholds & selection
      if (!is.null(safe_val("selected_lipid"))) {
        updateSelectizeInput(session, "selected_lipid",
                             selected = safe_val("selected_lipid"))
      }
      if (!is.null(safe_val("logFC_input")))
        updateNumericInput(session, "logFC_input",
                           value = safe_val("logFC_input"))
      if (!is.null(safe_val("p_value_max")))
        updateNumericInput(session, "p_value_max",
                           value = safe_val("p_value_max"))
      if (!is.null(safe_val("p_value_adj")))
        updateNumericInput(session, "p_value_adj",
                           value = safe_val("p_value_adj"))
      if (!is.null(safe_val("min_lipids_per_class")))
        updateNumericInput(session, "min_lipids_per_class",
                           value = safe_val("min_lipids_per_class"))
      
      # logFC scale & display
      if (!is.null(safe_val("selected_logfc_sclae_bar")))
        updateRadioButtons(session, "selected_logfc_sclae_bar",
                           selected = safe_val("selected_logfc_sclae_bar"))
      if (!is.null(safe_val("logFC_scale_manual")))
        updateNumericInput(session, "logFC_scale_manual",
                           value = safe_val("logFC_scale_manual"))
      if (!is.null(safe_val("split_screen")))
        updateCheckboxInput(session, "split_screen",
                            value = safe_val("split_screen"))
      if (!is.null(safe_val("show_grid")))
        updateCheckboxInput(session, "show_grid",
                            value = safe_val("show_grid"))
      
      # significance overlay
      if (!is.null(safe_val("outline_sig")))
        updateCheckboxInput(session, "outline_sig",
                            value = safe_val("outline_sig"))
      if (!is.null(safe_val("sig_shape")))
        updateRadioButtons(session, "sig_shape",
                           selected = safe_val("sig_shape"))
      if (!is.null(safe_val("sig_metric")))
        updateRadioButtons(session, "sig_metric",
                           selected = safe_val("sig_metric"))
      if (!is.null(safe_val("sig_threshold")))
        updateNumericInput(session, "sig_threshold",
                           value = safe_val("sig_threshold"))
      if (!is.null(safe_val("sig_outline_color")))
        colourpicker::updateColourInput(session, "sig_outline_color",
                                        value = safe_val("sig_outline_color"))
      if (!is.null(safe_val("sig_outline_size")))
        updateNumericInput(session, "sig_outline_size",
                           value = safe_val("sig_outline_size"))
      
      # color palette & tiles
      if (!is.null(safe_val("palette_preset")))
        updateSelectInput(session, "palette_preset",
                          selected = safe_val("palette_preset"))
      if (!is.null(safe_val("low_color")))
        colourpicker::updateColourInput(session, "low_color",
                                        value = safe_val("low_color"))
      if (!is.null(safe_val("mid_color")))
        colourpicker::updateColourInput(session, "mid_color",
                                        value = safe_val("mid_color"))
      if (!is.null(safe_val("high_color")))
        colourpicker::updateColourInput(session, "high_color",
                                        value = safe_val("high_color"))
      if (!is.null(safe_val("panel_bg_color")))
        colourpicker::updateColourInput(session, "panel_bg_color",
                                        value = safe_val("panel_bg_color"))
      if (!is.null(safe_val("strip_bg_color")))
        colourpicker::updateColourInput(session, "strip_bg_color",
                                        value = safe_val("strip_bg_color"))
      if (!is.null(safe_val("strip_text_color")))
        colourpicker::updateColourInput(session, "strip_text_color",
                                        value = safe_val("strip_text_color"))
      if (!is.null(safe_val("tile_border_show")))
        updateCheckboxInput(session, "tile_border_show",
                            value = safe_val("tile_border_show"))
      if (!is.null(safe_val("tile_border_color")))
        colourpicker::updateColourInput(session, "tile_border_color",
                                        value = safe_val("tile_border_color"))
      if (!is.null(safe_val("tile_border_size")))
        updateNumericInput(session, "tile_border_size",
                           value = safe_val("tile_border_size"))
      
      # text sizes & fonts
      if (!is.null(safe_val("strip_text_size")))
        updateNumericInput(session, "strip_text_size",
                           value = safe_val("strip_text_size"))
      if (!is.null(safe_val("axis_text_x_size")))
        updateNumericInput(session, "axis_text_x_size",
                           value = safe_val("axis_text_x_size"))
      if (!is.null(safe_val("axis_text_y_size")))
        updateNumericInput(session, "axis_text_y_size",
                           value = safe_val("axis_text_y_size"))
      if (!is.null(safe_val("axis_title_size")))
        updateNumericInput(session, "axis_title_size",
                           value = safe_val("axis_title_size"))
      if (!is.null(safe_val("legend_title_size")))
        updateNumericInput(session, "legend_title_size",
                           value = safe_val("legend_title_size"))
      if (!is.null(safe_val("legend_text_size")))
        updateNumericInput(session, "legend_text_size",
                           value = safe_val("legend_text_size"))
      if (!is.null(safe_val("axis_title_x_size")))
        updateNumericInput(session, "axis_title_x_size",
                           value = safe_val("axis_title_x_size"))
      if (!is.null(safe_val("axis_title_y_size")))
        updateNumericInput(session, "axis_title_y_size",
                           value = safe_val("axis_title_y_size"))
      if (!is.null(safe_val("plot_title_size")))
        updateNumericInput(session, "plot_title_size",
                           value = safe_val("plot_title_size"))
      if (!is.null(safe_val("plot_subtitle_size")))
        updateNumericInput(session, "plot_subtitle_size",
                           value = safe_val("plot_subtitle_size"))
      if (!is.null(safe_val("plot_subtitle_bold")))
        updateCheckboxInput(session, "plot_subtitle_bold",
                            value = safe_val("plot_subtitle_bold"))
      if (!is.null(safe_val("font_family")))
        updateSelectInput(session, "font_family",
                          selected = safe_val("font_family"))
      
      # axis / legend
      if (!is.null(safe_val("axis_text_x_angle")))
        updateNumericInput(session, "axis_text_x_angle",
                           value = safe_val("axis_text_x_angle"))
      if (!is.null(safe_val("barwidth")))
        updateNumericInput(session, "barwidth",
                           value = safe_val("barwidth"))
      if (!is.null(safe_val("barheight")))
        updateNumericInput(session, "barheight",
                           value = safe_val("barheight"))
      if (!is.null(safe_val("legend_side")))
        updateSelectInput(session, "legend_side",
                          selected = safe_val("legend_side"))
      
      # labels & titles
      if (!is.null(safe_val("x_axis_label")))
        updateTextInput(session, "x_axis_label",
                        value = safe_val("x_axis_label"))
      if (!is.null(safe_val("y_axis_label")))
        updateTextInput(session, "y_axis_label",
                        value = safe_val("y_axis_label"))
      if (!is.null(safe_val("main_title")))
        updateTextInput(session, "main_title",
                        value = safe_val("main_title"))
      if (!is.null(safe_val("sub_title")))
        updateTextInput(session, "sub_title",
                        value = safe_val("sub_title"))
      if (!is.null(safe_val("group_join")))
        updateTextInput(session, "group_join",
                        value = safe_val("group_join"))
      if (!is.null(safe_val("auto_subtitle_groups")))
        updateCheckboxInput(session, "auto_subtitle_groups",
                            value = safe_val("auto_subtitle_groups"))
      if (!is.null(safe_val("title_hjust")))
        updateNumericInput(session, "title_hjust",
                           value = safe_val("title_hjust"))
      if (!is.null(safe_val("subtitle_hjust")))
        updateNumericInput(session, "subtitle_hjust",
                           value = safe_val("subtitle_hjust"))
      
      # layout
      if (!is.null(safe_val("facet_cols")))
        updateNumericInput(session, "facet_cols",
                           value = safe_val("facet_cols"))
      if (!is.null(safe_val("plot_width_px")))
        updateNumericInput(session, "plot_width_px",
                           value = safe_val("plot_width_px"))
      if (!is.null(safe_val("plot_height_px")))
        updateNumericInput(session, "plot_height_px",
                           value = safe_val("plot_height_px"))
      if (!is.null(safe_val("panel_spacing_mm")))
        updateNumericInput(session, "panel_spacing_mm",
                           value = safe_val("panel_spacing_mm"))
      if (!is.null(safe_val("margin_top_mm")))
        updateNumericInput(session, "margin_top_mm",
                           value = safe_val("margin_top_mm"))
      if (!is.null(safe_val("margin_right_mm")))
        updateNumericInput(session, "margin_right_mm",
                           value = safe_val("margin_right_mm"))
      if (!is.null(safe_val("margin_bottom_mm")))
        updateNumericInput(session, "margin_bottom_mm",
                           value = safe_val("margin_bottom_mm"))
      if (!is.null(safe_val("margin_left_mm")))
        updateNumericInput(session, "margin_left_mm",
                           value = safe_val("margin_left_mm"))
    })
    
    
    
    
    # Message shown when hovering over Original data and merged data.
    observe({
      addTooltip(session, "selected_dataset", 
                 "Choose 'Original Data' to work with the data as it was initially collected. Select 'Merged Data' for a combined dataset.", 
                 placement = "bottem", 
                 trigger = "hover")
    })
    
    # Message shown when hovering over logFC. 
    observe({
      addTooltip(session, "logFC_input_ui", 
                 "Displays lipids where the absolute value of logFC is greater than or equal to the threshold. For example, entering '1' will include lipids with logFC ≥ 1 or logFC ≤ -1.", 
                 placement = "bottom", 
                 trigger = "hover")
    })
    
    
    # Tooltip for p-value max input
    observe({
      addTooltip(session, "p_value_max_ui", 
                 "Displays lipids where the p-value is less than or equal to the threshold.", 
                 placement = "bottom", 
                 trigger = "hover")
    })
    
    # Tooltip for p-value adjusted input
    observe({
      addTooltip(session, "p_value_adj", 
                 "Displays lipids where the adjusted p-value is less than or equal to the threshold.", 
                 placement = "bottom", 
                 trigger = "hover")
    })
    
    # Tooltip for minimum lipids per class input
    observe({
      addTooltip(session, "min_lipids_per_class_ui", 
                 "Includes lipid classes that have at least the specified minimum number of lipids.", 
                 placement = "bottom", 
                 trigger = "hover")
    })
    
    # Tooltip for logFC scale bar input
    observe({
      addTooltip(session, "selected_logfc_sclae_bar", 
                 "Select 'Dynamic' to use the dynamic range of logFC values, 
                   min and max of the logFC will then be set as the scale bar. 
                   Select 'Manual' to specify a custom range.", 
                 placement = "bottom", 
                 trigger = "hover")
    })
    
    # Shows which groups is selected for the logFC calculation in the 'Heatmap visualization' tab.
    output$selected_groups_heatmap_text <- renderUI({
      req(input$selected_group_for_numerator, input$selected_group_for_denominator)
      withMathJax(HTML(paste0(
        "<p>Data being compared is:</p>",
        "$$\\log_2\\left( \\frac{\\text{mean of group ", input$selected_group_for_numerator,
        "} + 10^{-6}}{\\text{mean of group ", input$selected_group_for_denominator,
        "} + 10^{-6}} \\right)$$"
      )))
    })
    
    
    
    
    
  })
}) # This finishes the first 'observeEvent' when 'Run data processing' is clicked


# Generate a nice summary table for the modal
output$cleaning_summary_ui <- renderUI({
  req(values$cleaning_stats)
  stats <- values$cleaning_stats
  
  # Calculate percentages
  pct_invalid <- round((stats$invalid / stats$start) * 100, 1)
  pct_final   <- round((stats$final / stats$start) * 100, 1)
  
  tagList(
    h4("Data Processing Statistics"),
    HTML(paste0(
      "<table class='table table-condensed' style='width:100%; font-size:14px;'>",
      "<tr style='background-color: #f9f9f9;'><td><strong>Total Imported Rows</strong></td><td style='text-align:right;'><strong>", stats$start, "</strong></td></tr>",
      
      "<tr><td>Rows Removed (Invalid Format)</td><td style='text-align:right; color:red;'>-", stats$invalid, " (", pct_invalid, "%)</td></tr>",
      
      "<tr><td><i>Rows Remaining after Cleaning</i></td><td style='text-align:right;'>", stats$valid, "</td></tr>",
      
      # Only show merge loss if it actually happened
      if(stats$merged_loss > 0) {
        paste0("<tr><td>Rows Collapsed (Isobaric Merging)</td><td style='text-align:right; color:orange;'>-", stats$merged_loss, "</td></tr>")
      } else { "" },
      
      "<tr style='border-top: 2px solid #333; background-color: #e6eef7;'><td><strong>Final Rows Available</strong></td><td style='text-align:right;'><strong>", stats$final, " (", pct_final, "%)</strong></td></tr>",
      "</table>"
    ))
  )
})

# Outside of the observeEvent, based on whether runProcessClicked is TRUE or FALSE, the message display will be placed on this: 
# For the first message, which is placed in the 'Lipid Heatmap' tab.
output$table_message_1 <- renderUI({
  if (!values$runProcessClicked) {
    HTML('<p>Make sure sequences file is uploaded, when uploaded: Press "Run Data Processing" to get a display of data</p>')
  }
})

# For the second message, which is placed in the 'Table' tab.
output$table_message_2 <- renderUI({
  if (!values$runProcessClicked) {
    HTML('<p>Make sure sequences file is uploaded, when uploaded: Press "Run Data Processing" to get a display of data</p>')
  }
})

# Outside of the observeEvent, so the message both are shown before and after runProcessClicked is clicked. 
observe({
  addTooltip(session, "selected_dataset", 
             "Choose 'Original Data' to work with the data as it was initially collected, all sum isoforms is within the data. Select 'Merged Data' for a combined dataset, meaning isoforms lipids are summed togehter.", 
             placement = "bottom", 
             trigger = "hover")
})

# User guide inside 'Heatmap'
observeEvent(input$show_lipid_info, {
  showModal(modalDialog(
    title = "Lipid Summary",
    textOutput("lipid_total"),  # Display the total number of lipids
    
    
    dataTableOutput("lipid_group_count"),  # Lipid group count table
    
    easyClose = TRUE,
    footer = modalButton("Close")
  ))
})

observeEvent(input$show_lipid_cal, {
  showModal(modalDialog(
    title = "Lipid Calculation",
    div(
      tags$b("logFC: "), 
      withMathJax("$$\\log_2\\left( \\frac{\\text{mean of group numerator} + 10^{-6}}{\\text{mean of group denominator} + 10^{-6}} \\right)$$")
    ),
    div(
      tags$b("logFC Explanation: "), 
      "+10^-6 is added to avoid division by zero. The logFC is calculated for each lipid, 
      comparing the mean of the numerator group to the mean of the denominator group, 
      ensuring all sample values are taken into account."
    ),
    div(
      tags$b("p-value: "), 
      "The p-value is calculated using Welch’s t-test, comparing the raw data between the two groups for each lipid. 
      The data used is not filtered by p-value or logFC thresholds beforehand, ensuring an unbiased comparison."
    ),
    div(
      tags$b("Adjusted p-value (p-adj): "), 
      "P-values are adjusted for multiple comparisons using the Benjamini-Hochberg (BH) method. 
      This controls the false discovery rate (FDR), reducing the likelihood of false positives when testing multiple lipids simultaneously."
    ),
    
    div(
      tags$b("Packages Used: "),
      "This analysis utilizes the 'lipidomeR' package version 0.1.2."
    ),
    easyClose = TRUE,
    footer = modalButton("Close")
  ))
})

# Render the datatable of removed_data
output$lipid_remove_table <- DT::renderDataTable({
  req(input$show_lipid_remove, values$removed_data) # Ensures button has been clicked and data exists
  DT::datatable(
    values$removed_data,
    options = list(pageLength = 5, autoWidth = TRUE, scrollX = TRUE)
  )
})

# Show the "Wise" modal dialog
observeEvent(input$show_lipid_remove, {
  showModal(modalDialog(
    title = "Lipid Filtration Summary", # Fixed spelling
    size = "l",                         # Made it 'large' to fit content better
    
    # 1. The new wise summary table
    uiOutput("cleaning_summary_ui"),
    
    hr(),
    
    # 2. The details of what was removed
    h4("Details: Rows Removed (Invalid Format)"),
    p("The following rows were removed because their names in the first column did not match the required 'Class(C:D)' format (e.g., 'PC(34:1)')."),
    DT::dataTableOutput("lipid_remove_table"),
    
    easyClose = TRUE,
    footer = modalButton("Close")
  ))
})

observeEvent(input$df_explain, {
  showModal(modalDialog(
    title = "Lipid Data Frames",
    div(
      tags$b("Original vs Merged Data"),
      tags$p("In 'Original data', every molecular feature is shown separately. This includes isobaric species—molecules that share the same mass but differ in structure. Lipidomics datasets often contain thousands of lipids, where even subtle differences (such as double bond positions) matter."),
      tags$p("For example, consider the lipids:"),
      tags$ul(
        tags$li("DG(18:1(9Z)/20:3(8Z,11Z,14Z)/0:0)[iso2]"),
        tags$li("DG(18:0/20:4(5Z,8Z,11Z,14Z)/0:0)[iso2]"),
        tags$p("These lipids share the same sum composition (DG(38:4)) but differ in structure and abundance.")
      ),
      
      tags$p("Although both have the same mass (m/z 644.537975), their structures and abundances differ. In 'Original data', these would appear as DG_1(38:4) and DG_2(38:4), maintaining each unique identity since the program only supports the X(C:D) notation."),
      
      tags$p("In contrast,'Merged data'combines these isobaric species into a single entry (DG(38:4)), summing their abundances and removing the structural distinctions that existed in the original data."),
      
      tags$p("References:"),
      tags$ul(
        tags$li(tags$a(href = "https://www.lipidmaps.org/databases/lmsd/LMGL02010110", "DG 18:1_20:3")),
        tags$li(tags$a(href = "https://www.lipidmaps.org/databases/lmsd/LMGL02010111", "DG 18:0_20:4"))
      ),
      tags$p("To learn more about these formats, please visit the ## WILL BE ADDED I MAIN_SERVER",
             tags$a(href = "http://documentation_link_here_(will_come_in_main_server", "documentation page"), "."),
    ),
    easyClose = TRUE,
    footer = modalButton("Close")
  ))
})

observeEvent(input$lipid_contact, {
  showModal(
    modalDialog(
      title = "Contact Information",
      div(
        tags$p(tags$b("Lipid Heatmap")),
        tags$p("This program is written in R and Shiny, and is designed to help researchers and students to analyze lipidomics data."),
        tags$p("It is a part of the Bachelor's thesis project at the University of Southern Denmark, and is developed by Jakob R. Hartvig. BS.c. in Biomedicine"), 
        tags$p("Multiple datasets have been tested—both newly generated and previously published—to ensure that the program works as intended."),
        
        tags$p("Feel free to use, modify, and distribute this tool under the terms of its open-source license. If you publish results obtained with this tool, 
                 a citation or acknowledgment would be greatly appreciated."), 
        tags$p("The program is open-source and can be found on Metabolink GitHub:"),
        tags$p(tags$b("If you have any questions, bug finds or concerns, please get in touch:")),
        tags$ul(
          tags$li("Jakob R. Hartvig"),
          tags$li("Email: ", tags$a(href = "mailto:jahar21@student.sdu.dk", "jahar21@student.sdu.dk")),
          tags$li("If there is no response within two to three business days, please try: ",
                  tags$a(href = "mailto:jakobhartvigprivat@gmail.com", "jakobhartvigprivat@gmail.com"),
          ),
          tags$p("I welcome suggestions, feedback, and contributions from the community!"),
          tags$p("Thank you for using Lipid Heatmap!"),
          tags$p("Best regards, Jakob R. Hartvig"),
          
          
          
          tags$svg(
            width = "120px", height = "140px", viewBox = "0 0 120 140",
            style = "display:block; margin:auto;",
            
            # Define a linear gradient for the droplet (from bright yellow at the top to warm pink at the bottom)
            tags$defs(
              tags$linearGradient(
                id = "funGrad", x1 = "0%", y1 = "0%", x2 = "0%", y2 = "100%",
                tags$stop(offset = "0%", style = "stop-color:#fdf21c;stop-opacity:1"), # Bright yellow top
                tags$stop(offset = "100%", style = "stop-color:#ff7d7d;stop-opacity:1") # Warm pink bottom
              )
            ),
            
            # Draw a droplet shape
            tags$path(
              d = "M60,10 
         C 20,10 0,50 60,120 
         C120,50 100,10 60,10 Z",
              fill = "url(#funGrad)",
              stroke = "#333",
              `stroke-width` = 2
            ),
            
            # Add two eyes (small black circles)
            tags$circle(cx="50", cy="50", r="4", fill="#333"),
            tags$circle(cx="70", cy="50", r="4", fill="#333"),
            
            # Add a smiling mouth using a path (a simple arc)
            tags$path(
              d = "M45,65 Q60,80 75,65", 
              fill = "transparent",
              stroke = "#333",
              `stroke-width` = 3,
              `stroke-linecap` = "round"
            ),
            
            # Add the text "Lipid" below the droplet in a bright color
            tags$text(
              "LIPIDS",
              x = "50%", y = "131",
              `text-anchor` = "middle",
              fill = "#ff4f4f",
              style = "font-family:sans-serif; font-size:22px; font-weight:bold;"
            )
          )
          
          
          
        ),
      ),
      easyClose = TRUE,
      footer = modalButton("Close")
    )
  )
}) #### End of all lipid heatmap code for server.r 
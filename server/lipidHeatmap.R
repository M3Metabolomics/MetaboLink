  # Values is used for message display before and after data process
  values <- reactiveValues(runProcessClicked = FALSE)
  
  # When bottom clicked in interface, all the following will be processed
  observeEvent(input$run_process, {
    values$runProcessClicked <- TRUE
    
    # Required data files are loaded
    sequence <- rv$sequence[[rv$activeFile]]
    data <- rv$data[[rv$activeFile]]
    
    # Capture the number of rows before filtering, used to compare with the number of rows after filtering
    number_of_rows_before <- nrow(data)
    
    
    ###############
    # Data cleaning
    ###############
    
    # Removes noise, keeping only the name and length (e.g., "CAR 14:1'CAR'[M+H]+" becomes "CAR 14:1")
    data[, 1] <- sapply(data[, 1], extract_pattern)
    
    # Apply the `reposition_ether_lipids` function to change the posistion of ether lipids (O- or P- prefixed lipids). e.g. LP(O-18:1) --> LP-O(18:1)
    data[, 1] <- sapply(data[, 1], reposition_ether_lipids)
    
    # Apply the `clean_lipid_name` function to remove content after semicolons inside parentheses
    data[, 1] <- sapply(data[, 1], clean_lipid_name)
    
    # Apply the `format_strings` function to ensure the lipid names are properly formatted, adding parentheses around numbers if necessary
    data[, 1] <- sapply(data[, 1], format_strings)
    
    # Used to show which data is being filtraed away from the data frame
    removed_data <<- remove_patterned_rows(data)
    
    # Function to filter rows based on the specified pattern, meaning removes any data that are not on X(C:D) format.
    data <- filter_data_by_pattern(data)
    
    # Before the if-block, initialize merged_data_info
    merged_data_info <- NULL
    
    
    # This will make it possible to switch between original data and merged data. OG data: using _1 _2 ... _n. Merged will sum the values of the "duplicated" data. 
    if (input$selected_dataset == "original") {
      data <- unique_compound_names(data)
      
      # Merge duplicates
    }  else if (input$selected_dataset == "merged") {
      
      # Creating info to display inside 'Table of Heatmap' about which lipids are merged. 
      merged_data_info <- merged_info_function(data)
      
      # Removes anything that are not part of the data of the samples and name. So it it possible to sum the rows. 
      data <- data[, sequence[ , 'labels'] %in% c("Name","Sample")]
      sequence <- sequence[sequence[ , 'labels'] %in% c("Name","Sample"), ]
      
      # Sums the duplicates (isoforms of the lipids)
      data <- merge_duplicates(data)
    }
    
    
    # Capture the number of rows after filtering
    number_of_rows_after <- nrow(data)
    
    # Calculate the number of rows removed
    rows_removed <- number_of_rows_before - number_of_rows_after
    output$rows_removed_text <- renderText({
      paste("Rows removed after data cleaning are:", rows_removed, ". The removal is be due to the names in the first column of the data file not being in the X(C:D) format. 
            Keep in mind, that the merged data will also count as a removed row, but they will not show up in the table down below.")
    })
    
    
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
      grouped_data_frames <<- create_grouped_data_frames(sequence, data)
      
      compound_names <- data[[1]]  # Extract the first column which contains compound names
      
      # Assuming that each grouped data frame has rows in the same order as "data"
      for (i in seq_along(grouped_data_frames)) {
        grouped_data_frames[[i]] <- cbind(Compound_Name = compound_names, grouped_data_frames[[i]])
      }
      
      # Extract unique group names from sequence
      unique_group_names <- unique(sequence[sequence$labels == "Sample", "group"])
      
      # Check if lengths of group names and grouped data frames match
      if (length(unique_group_names) == length(grouped_data_frames)) {
        # Apply the actual group names from the sequence file
        names(grouped_data_frames) <- unique_group_names
      } else {
        # If there's a mismatch, fallback to naming with numbers as before
        names(grouped_data_frames) <- paste("Group", seq_along(grouped_data_frames))
      }
      
      # Render the UI for group selection.
      # Render the UI for group selection.
      output$select_group_ui_heatmap <- renderUI({
        column(
          title = "Select groups for Heatmap test",
          width = 12,
          tagList(
            selectInput("selected_group_for_numerator", "Select Group for numerator:",
                        choices = names(grouped_data_frames),
                        selected = names(grouped_data_frames)[1] # Default to the first group
            ),
            selectInput("selected_group_for_denominator", "Select Group for denominator:",
                        choices = names(grouped_data_frames),
                        selected = names(grouped_data_frames)[2] # Default to the first group
            )
          )
        )
      })
      
      
      
      
      # Create interactive table for selected numerator group. Displayed at the starting page of the 'Lipid Heatmap'
      output$numerator_group_table <- DT::renderDataTable({
        req(input$selected_group_for_numerator) 
        # Create a copy of the data for display purposes
        display_data <- grouped_data_frames[[input$selected_group_for_numerator]]
        
        DT::datatable(
          display_data,  
          options = list(scrollX = TRUE)  # Enable horizontal scrolling
        )
      })
      
      # Create interactive table for selected denominator group
      output$denominator_group_table <- DT::renderDataTable({
        req(input$selected_group_for_denominator)  
        display_data <- grouped_data_frames[[input$selected_group_for_denominator]]
        
        DT::datatable(
          display_data,  
          options = list(scrollX = TRUE)  
        )
      })
      
      
      # Interface of selections of lipids to display
      output$select_lipid_ui <- renderUI({
        # Extract the lipid names from first column of the file 'data'
        lipid_names <<- group_lipids_by_group(data)
        
        selectizeInput("selected_lipid", "Select lipid(s) to display:",
                       choices = c("All", unique(lipid_names$group)),
                       multiple = TRUE,
                       options = list(placeholder = 'Choose lipids...',
                                      onInitialize = I('function() { this.setValue(""); }')))
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
        
        numerator_data <- grouped_data_frames[[input$selected_group_for_numerator]]
        denominator_data <- grouped_data_frames[[input$selected_group_for_denominator]]
        
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
            t_test_result <- t.test(num_values, denom_values)
            p_values[i] <- t_test_result$p.value
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
        numerator_data <- grouped_data_frames[[input$selected_group_for_numerator]]
        denominator_data <- grouped_data_frames[[input$selected_group_for_denominator]]
        
        # Removes the fir column of the data, as it is not needed for the calculation.
        numerator_data <- numerator_data[, -1]
        denominator_data <- denominator_data[, -1]
        
        numerator_data_means <- rowMeans(numerator_data, na.rm = TRUE) # Makes sure that the calculation is done even NA values are present. # "drop = FALSE" was removed to avoid the error message.
        denominator_data_means <- rowMeans(denominator_data, na.rm = TRUE) # "drop = FALSE" was removed to avoid the error message.
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
          logFC <- logFC[lipid_names$group %in% input$selected_lipid, ]
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
      
      # Reactive expression to filter data based on p-value and logFC thresholds, plus the amount of lipids within their class.
      reactiveFilteredData <- reactive({
        logFC_data <- reactiveLogFC()  
        p_value_data <- reactiveP_value()
        
        filtered_data <- merge(logFC_data, p_value_data, by = "Compound_Name")
        
        # Apply initial filtering criteria
        filtered_data <- filtered_data %>%
          filter(!is.na(p_value) & p_value <= input$p_value_max) %>%
          filter(abs(logFC) >= input$logFC_input)
        
        # Check if filtered_data has rows after initial filtering
        if (nrow(filtered_data) == 0) {
          # No data to proceed with, return empty data frame
          return(filtered_data)
        }
        
        # Get lipid class mapping for all lipids in the original data
        names.mapping.all <- map_lipid_names(x = filtered_data[[1]])
        
        # Compute counts per class from the original data
        class_counts <- table(names.mapping.all$Class)
        class_counts_df <- as.data.frame(class_counts)
        colnames(class_counts_df) <- c("Class", "Count")
        
        # Set the threshold from user input
        min_lipids_per_class <- input$min_lipids_per_class
        
        # Get classes that meet the threshold
        classes_to_keep <- class_counts_df$Class[class_counts_df$Count >= min_lipids_per_class]
        
        # Map the class information to the filtered_data
        # Create a data frame of lipid names and class
        lipid_class_df <- data.frame(Compound_Name = filtered_data[[1]], Class = names.mapping.all$Class)
        
        # Merge the class information with filtered_data
        filtered_data <- merge(filtered_data, lipid_class_df, by = "Compound_Name", all.x = TRUE)
        
        # Now filter based on classes_to_keep
        filtered_data <- filtered_data[filtered_data$Class %in% classes_to_keep, ]
        
        # Check again if filtered_data has rows after class filtering
        if (nrow(filtered_data) == 0) {
          # No data to proceed with, return empty data frame
          return(filtered_data)
        }
        
        return(filtered_data)
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
      
      
      # Heatmap plot loading
      
      # Reactive expression for heatmap plot data and height
      heatmap_plot_data <- reactive({
        # Used the data frame which is adjusted to the user input.
        filtered_data <- reactiveFilteredData()
        
        # Ensure the data is not NULL and has rows to plot
        req(nrow(filtered_data) > 0)
        
        # Map the lipid names (make sure it returns all necessary columns, including Lipid_Class)
        # This is a function from lipidomeR package. 
        names.mapping <- map_lipid_names(x = filtered_data$Compound_Name)
        
        # Calculate the number of unique classes
        num_classes <- length(unique(names.mapping$Class))
        
        # Calculate number of rows in the facets
        ncol_facets <- 3  # Adjust based on your facet_wrap setting
        num_rows <- ceiling(num_classes / ncol_facets)
        
        # Set base height per row
        height_per_row <- 300  # Adjust as needed to match bubble plot scaling
        
        # Calculate total plot height
        total_plot_height <- num_rows * height_per_row
        
        # Ensure minimum and maximum height limits
        total_plot_height <- max(total_plot_height, 600)   # Minimum height
        total_plot_height <- base::min(total_plot_height, 4000)  # Maximum height
        
        # Calculate logFC range dynamically
        logFC_range <- range(filtered_data$logFC, na.rm = TRUE)
        
        # Adjust fill.limits based on user input
        if (input$selected_logfc_sclae_bar == "Manual") {
          manual_limit <- input$logFC_scale_manual
          fill.limits <- c(-manual_limit, manual_limit)
        } else {
          # Use dynamic range
          fill.limits <- c(base::min(logFC_range), max(logFC_range))
        }
        
        heatmap_plot <- heatmap_lipidome(
          x = filtered_data[ , c("Compound_Name", "logFC")],
          names.mapping = names.mapping,
          class.facet = "wrap",
          x.names = "Compound_Name",
          fill.limits = fill.limits,  
          melt.value.name = "logFC",
          scales = "free"
        ) +
          scale_fill_gradient2(
            low = input$low_color,
            mid = input$mid_color,
            high = input$high_color,
            limits = fill.limits,  
            space = "Lab",
            name = "logFC",
            guide = guide_colorbar(
              barwidth = input$barwidth, # Adjust width of color bar
              barheight = input$barheight  # Adjust height of color bar
            )
          ) +
          facet_wrap(~ Class, scales = "free", ncol = ncol_facets) +
          theme(
            panel.background = element_rect(fill = input$panel_bg_color, color = "white"),
            strip.background = element_rect(fill = input$strip_bg_color, color = "white"),
            strip.text = element_text(color = input$strip_text_color, face = "bold", size = input$strip_text_size),
            axis.text.x = element_text(
              angle = input$axis_text_x_angle, 
              size = input$axis_text_x_size
            ),
            axis.text.y = element_text(size = input$axis_text_y_size),
            axis.title = element_text(size = input$axis_title_size),
            legend.title = element_text(size = input$legend_title_size),
            legend.text = element_text(size = input$legend_text_size),
            axis.title.x = element_text(size = input$axis_title_x_size),
            axis.title.y = element_text(size = input$axis_title_y_size),
            plot.title = element_text(size = input$plot_title_size, face = "bold")
          ) +
          labs(
            x = input$x_axis_label,
            y = input$y_axis_label,
            title = input$main_title, 
          ) 
        
        
        # Conditionally remove grid lines based on user input
        if (!input$show_grid) {
          heatmap_plot <- heatmap_plot +
            theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
            )
        }
        
        # Return the plot and height
        list(plot = heatmap_plot, height = total_plot_height)
      })
      
      # Render the heatmap plot
      output$heatmapPlot <- renderPlot({
        heatmap_plot_data()$plot
      }, height = function() {
        heatmap_plot_data()$height
      })
      
      
      
      # IF user want to split the screen, the following will be shown. Heatmap and table side by side.
      output$visualization_ui <- renderUI({
        if (input$split_screen) {
          # Show both heatmap and table side by side
          fluidRow(
            column(width = 6,
                   div(style = "width: 100%; height: 800px; overflow-y: scroll;",  
                       withSpinner(plotOutput("heatmapPlot", width = "100%", height = "800px"))
                   )
            ),
            column(width = 6,
                   div(style = "width: 100%; height: 800px; overflow-y: scroll;",  
                       withSpinner(dataTableOutput("pValueTable_2"))
                   )
            )
          )
        } else {
          # Show only heatmap
          div(style = "width: 100%; height: 2000px; overflow-y: scroll;",  
              withSpinner(plotOutput("heatmapPlot", width = "100%", height = "2000px"))
          )
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
        
        # If 'merged_data_info' is required only in "merged" scenario, don't use req(merged_data_info)
        # Instead, just check if it's not NULL before merging later.
        
        # Ensure 'Compound_Name' exists in 'data'
        if (!"Compound_Name" %in% colnames(data)) {
          data$Compound_Name <- data[, 1]  # Assuming the first column contains compound names
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
      
      
      
      # Observe the action button to open the modal dialog
      observeEvent(input$download_heatmap_btn, {
        showModal(modalDialog(
          title = "Download Heatmap Image",
          selectInput("modal_image_format", "Select Image Format", choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
          numericInput("modal_image_dpi", "Image Resolution (DPI)", value = 300, min = 72, step = 72),
          footer = tagList(
            modalButton("Cancel"),
            downloadButton("modal_download_heatmap", "Download")
          )
        ))
      })
      
      
      
      
      # Reactive expression to calculate plot dimensions
      plot_dimensions <- reactive({
        # Take the input from user in interface and change p-value and logFC
        filtered_data <- reactiveFilteredData()
        
        # Ensure the data is not NULL and has rows to plot
        req(nrow(filtered_data) > 0)
        
        # Map the lipid names (make sure it returns all necessary columns, including Lipid_Class)
        names.mapping <- map_lipid_names(x = filtered_data$Compound_Name)
        
        # Calculate the number of unique classes
        num_classes <- length(unique(names.mapping$Class))
        
        # Calculate number of rows in the facets
        ncol_facets <- 3  # Adjust based on your facet_wrap setting
        num_rows <- ceiling(num_classes / ncol_facets)
        
        # Set base dimensions per row/column
        height_per_row <- 3  # in inches
        width_per_col <- 4   # in inches
        
        # Calculate total plot dimensions
        total_height <- num_rows * height_per_row
        total_width <- ncol_facets * width_per_col
        
        # Ensure minimum and maximum limits
        total_height <- max(total_height, 6)   # Minimum height in inches
        total_height <- base::min(total_height, 40)  # Maximum height in inches
        total_width <- max(total_width, 8)     # Minimum width in inches
        total_width <- base::min(total_width, 40)    # Maximum width in inches
        
        # Return the dimensions
        list(height = total_height, width = total_width)
      })
      
      
      
      # Download handler for the heatmap plot from the modal dialog
      output$modal_download_heatmap <- downloadHandler(
        filename = function() {
          paste("heatmap_plot_", Sys.Date(), ".", input$modal_image_format, sep = "")
        },
        content = function(file) {
          # Generate the plot
          heatmap_plot <- heatmap_plot_data()$plot
          
          # Get plot dimensions
          dims <- plot_dimensions()
          total_height <- dims$height  # in inches
          total_width <- dims$width    # in inches
          
          # Save the plot with user-specified format and resolution
          ggsave(
            filename = file,
            plot = heatmap_plot,
            device = input$modal_image_format,
            width = total_width,    # Width in inches
            height = total_height,  # Height in inches
            units = "in",
            dpi = input$modal_image_dpi
          )
          
          # Close the modal dialog after download
          removeModal()
        }
      )
      
      
      
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
    req(input$show_lipid_remove) # Ensures button has been clicked
    DT::datatable(
      removed_data,
      options = list(pageLength = 5, autoWidth = TRUE, scrollX = TRUE)
    )
  })
  
  # Show the modal dialog containing the table when button is pressed
  observeEvent(input$show_lipid_remove, {
    showModal(modalDialog(
      title = "Lipid filtration summeray",
      textOutput("rows_removed_text"),
      dataTableOutput("lipid_remove_table"),
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
  })
  
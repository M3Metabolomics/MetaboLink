  observeEvent(input$run_circular_barplot, {
    # Ensure all required inputs are available.
    req(input$select_data_circular_barplot,
        input$top_x_cirbar,
        input$name_column_cirbar,
        input$group_column_cirbar,  
        input$group1_cirbar,
        input$group2_cirbar)
    
    if (!is.null(rv$activeFile)) {
      if (input$select_data_circular_barplot == "Unsaved data") {
        data <- rv$tmpData  # Use the temporary data
        seq <- rv$tmpSequence  # Use the temporary sequence
        dataset_name <- "Unsaved data"
      } else {
        # Get the index of the selected dataset
        sd <- which(rv$choices %in% input$select_data_circular_barplot)
        data <- rv$data[[sd]]  # Retrieve the selected dataset
        seq <- rv$sequence[[sd]]  # Retrieve the selected sequence
        dataset_name <- names(rv$data)[sd]  # Retrieve dataset name
      }
      
      # Retrieve input selections
      names_col   <- input$name_column_cirbar
      grouping    <- input$group_column_cirbar  # use the correct input id
      numerator      <- input$group1_cirbar
      denominator      <- input$group2_cirbar
      top_X <- input$top_x_cirbar
      feature <- input$feature_cirbar
      
      if (numerator == denominator) {
        sendSweetAlert(session, "Error", "Please select different groups for comparison.", type = "error")
        return()
      }
      
      seq_subset <- seq[seq[, "labels"] %in% c("Sample"), ]  # Restrict to "Sample" rows
      
      # if seq[,group] is numeric put group_ in front
      if (is.numeric(seq_subset[, "group"])) {
        seq_subset[, "group"] <- paste0("group_", seq_subset[, "group"])
        numerator <- paste0("group_", numerator)
        denominator <- paste0("group_", denominator)
      }
      
      # remove rows where data[, names_col] is empty or NA
      data <- data[!is.na(data[, names_col]) & data[, names_col] != "", ]
      
      data_subset <- data[, c(rownames(seq_subset)), drop = FALSE]  # Use row names of seq_subset to filter columns
      
      rownames(data_subset) <- make.unique(rep(data[, names_col], length.out = nrow(data_subset)))

      stat_results <- calculate_stats(data_subset, seq_subset, adjust.method = "BH")
      
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
      
      # assign the sub_df a name for debugging
      cirbar_df_name <- paste0(dataset_name, "_cirbar_df")
      assign(cirbar_df_name, sub_df)
      
      # Define custom colors
      custom_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                         "#D55E00", "#CC79A7", "#999999", "#8E44AD", "#1ABC9C",
                         "#2ECC71", "#3498DB", "#F39C12", "#E74C3C", "#95A5A6",
                         "#34495E", "#FB9A99", "brown", "#2980B9", "#C0392B")
      
      # column bind with dplyr 
      plotting_data <- bind_cols(
        data %>% select(all_of(names_col),all_of(grouping)),
        sub_df %>% select(FC, log2FC, p.value, p.adj, AveExpr, t, B)
      )
      
      # for each column return the min and max value in a 2xncol(data) matrix
      min_max_values <- sapply(plotting_data[,3:ncol(plotting_data)], function(x) c(base::min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
      rownames(min_max_values) <- c("min", "max")
      print(round(min_max_values, 2))
      
      plotting_data <- plotting_data %>%
        mutate(!!grouping := as.factor(.data[[grouping]])) %>%
        arrange(p.adj) %>% 
        # select top_X rows 
        slice_head(n = top_X) %>%
        arrange(.data[[grouping]],p.adj)
      
      # Add empty start and end rows
      # empty_start <- data.frame(matrix(0, nrow = 1, ncol = ncol(plotting_data)))
      # colnames(empty_start) <- colnames(plotting_data)
      # empty_start[1, 1:2] <- NA
      
      empty_end <- data.frame(matrix(0, nrow = 1, ncol = ncol(plotting_data)))
      colnames(empty_end) <- colnames(plotting_data)
      empty_end[1, 1:2] <- NA
      
      plotting_data <- rbind(plotting_data, empty_end)
      plotting_data <- plotting_data %>%
        mutate(id = row_number()) %>%
        relocate(id, .before = all_of(names_col)) # Reassign id after adding empty bars
      
      # Remove potential NA values in factor levels
      plotting_data[,grouping] <- factor(plotting_data[,grouping], exclude = NULL)
      
      # Prepare label data
      label_data <- plotting_data
      number_of_bar <- nrow(label_data)

      # Calculate angles and alignment
      label_data <- label_data %>%
        mutate(angle = 90 - 360 * (id - 0.5) / number_of_bar,
               hjust = ifelse(angle < -90, 1, 0),
               angle = ifelse(angle < -90, angle + 180, angle))
      
      valid_data <- plotting_data %>% filter(!is.na(.data[[feature]]))
      
      # Create 5 evenly spaced ticks from min to max
      tick_min <- base::min(plotting_data[[feature]], na.rm = TRUE)
      tick_max <- max(plotting_data[[feature]], na.rm = TRUE)
      ticks <- seq(tick_min, tick_max, length.out = 5)
      
      tick_1 <- ticks[1]  # same as min
      tick_2 <- ticks[2]  # 25% of the range
      tick_3 <- ticks[3]  # 50% of the range
      tick_4 <- ticks[4]  # 75% of the range
      tick_5 <- ticks[5]  # same as max
      
      # Generate the circular bar plot with correctly ordered
      cirbar_plot <- ggplot(plotting_data,
                            aes(x = as.factor(id),
                                y = .data[[feature]],
                                fill = .data[[grouping]])) +
        geom_bar(stat = "identity", alpha = 0.5) +
        
        # geom_segment(data = valid_data,
        #              aes(x = as.factor(id), xend = as.factor(id), y = tick_1, yend = tick_5),
        #              linetype = "dashed", color = "gray", linewidth = 0.3) +
        
        geom_segment(data = valid_data,
                     aes(x = base::min(id), xend = max(id), y = tick_1, yend = tick_1),
                     linetype = "dashed", color = "gray70", linewidth = 0.3) +
        geom_segment(data = valid_data,
                     aes(x = base::min(id), xend = max(id), y = tick_2, yend = tick_2),
                     linetype = "dashed", color = "gray70", linewidth = 0.3) +
        geom_segment(data = valid_data,
                     aes(x = base::min(id), xend = max(id), y = tick_3, yend = tick_3),
                     linetype = "dashed", color = "gray70", linewidth = 0.3) +
        geom_segment(data = valid_data,
                     aes(x = base::min(id), xend = max(id), y = tick_4, yend = tick_4),
                     linetype = "dashed", color = "gray70", linewidth = 0.3) +
        geom_segment(data = valid_data,
                     aes(x = base::min(id), xend = max(id), y = tick_5, yend = tick_5),
                     linetype = "dashed", color = "gray70", linewidth = 0.3) +
        
        annotate("text",
                 x = base::min(plotting_data$id),
                 y = c(tick_1, tick_2, tick_3, tick_4, tick_5),
                 label = c(as.character(round(tick_1, 2)),
                           as.character(round(tick_2, 2)),
                           as.character(round(tick_3, 2)),
                           as.character(round(tick_4, 2)),
                           as.character(round(tick_5, 2))),
                 color = "black", size = 5, angle = 0, fontface = "bold", hjust = 1) +
        
        ylim(base::min(valid_data[,feature]) - 0.25,
             max(valid_data[,feature])) +
        
        theme_minimal() +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.margin = unit(rep(0, 4), "cm")
        ) +
        coord_polar(start = 0) +
        scale_fill_manual(values = custom_colors, na.value = "gray90") +
        geom_text(data = label_data %>% filter(!is.na(.data[[names_col]])),
                  aes(x = as.factor(id),
                      y = tick_5,
                      label = .data[[names_col]],
                      hjust = hjust,
                      angle = angle),
                  color = "black",
                  fontface = "bold",
                  alpha = 1,
                  size = 3,
                  inherit.aes = FALSE)
      
      # Render the plotly plot
      output$circular_barplot <- renderPlot({
        cirbar_plot
      })
      
      # TODO: return plotting_data as DF to investigate the data used for the plot
      
    }
  })
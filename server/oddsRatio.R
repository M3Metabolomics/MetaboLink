  
  observeEvent({
    list(input$select_OR_data, input$OR_main_label)
  }, {
    req(input$select_OR_data)
    
    if (!is.null(rv$activeFile)) {
      # Retrieve the dataset and sequence based on selection
      if (input$select_OR_data == "Unsaved data") {
        data <- rv$tmpData    # Use temporary data
        seq  <- rv$tmpSequence
      } else {
        # Find the index of the selected dataset
        sd   <- which(rv$choices %in% input$select_OR_data)
        data <- rv$data[[sd]]
        seq  <- rv$sequence[[sd]]
      }
      
      # Extract column names from the selected dataset
      data_colnames <- colnames(data)
      
      # Define default values based on data columns
      main_default    <- if ("super_class" %in% data_colnames) "super_class" else if ("main_class" %in% data_colnames) "main_class" else data_colnames[1]
      sub_default     <- if ("sub_class" %in% data_colnames) "sub_class" else if ("Lipid.Abbreviation" %in% data_colnames) "Lipid.Abbreviation" else data_colnames[1]
      feature_default <- if ("Original.annotation" %in% data_colnames) "Original.annotation" else if ("Species.Name" %in% data_colnames) "Species.Name" else data_colnames[1]
      
      # Update the select inputs using the data columns;
      updateSelectInput(session, "OR_main_label",
                        choices  = data_colnames,
                        selected = if (input$OR_main_label %in% data_colnames) input$OR_main_label else main_default)
      
      updateSelectInput(session, "OR_sub_label",
                        choices  = data_colnames,
                        selected = if (input$OR_sub_label %in% data_colnames) input$OR_sub_label else sub_default)
      
      updateSelectInput(session, "OR_feature_label",
                        choices  = data_colnames,
                        selected = if (input$OR_feature_label %in% data_colnames) input$OR_feature_label else feature_default)
      
      # Determine which main label to use for computing groups:
      selected_main_label <- if (input$OR_main_label %in% data_colnames) {
        input$OR_main_label
      } else {
        main_default
      }
      
      # Extract unique groups from the determined main label column and remove empty strings
      groups <- sort(unique(data[[selected_main_label]]))
      groups <- groups[nzchar(groups)]
      
      # Create a choices vector with an "All" option as the first element.
      choices <- c("All", groups)
      
      # Update the selectize input for groups with "All" as the default selection.
      updateSelectizeInput(session, "selected_groups_OR",
                           choices  = choices,
                           selected = "All",
                           server   = TRUE)
    }
  })
  
  # OR plot
  observeEvent(input$run_OR_plot, {
    req(input$select_OR_data,
        input$group1_OR,
        input$group2_OR,
        input$OR_main_label,
        input$OR_sub_label,
        input$OR_feature_label,
        input$selected_groups_OR)
    
    if (!is.null(rv$activeFile)) {
      if (input$select_OR_data == "Unsaved data") {
        data <- rv$tmpData  
        seq  <- rv$tmpSequence     
        dataset_name <- "Unsaved data"
      } else {
        sd <- which(rv$choices %in% input$select_OR_data)
        data <- rv$data[[sd]]      
        seq  <- rv$sequence[[sd]]  
        dataset_name <- names(rv$data)[sd]
      }
      
      seq <- seq[seq$labels == "Sample", ]
      
      group1 <- input$group1_OR
      group2 <- input$group2_OR
      Main_label <- input$OR_main_label        
      Sub_label <- input$OR_sub_label          
      Feature_label <- input$OR_feature_label  
      Selected_Groups <- input$selected_groups_OR
      
      # When processing the selected groups:
      Selected_Groups <- if ("All" %in% input$selected_groups_OR) {
        data[[Main_label]] %>% unique() %>% sort()
      } else {
        input$selected_groups_OR
      }
      
      message("Inside Odds ratio")
      message("Group 1: ", group1)
      message("Group 2: ", group2)
      message("Main label: ", Main_label)
      message("Sub label: ", Sub_label)
      message("Feature label: ", Feature_label)
      message("Selected features: ", Selected_Groups)
      
      # Main_label, Sub_label and Feature_label must not be identical 
      if (Main_label == Sub_label || Main_label == Feature_label || Sub_label == Feature_label) {
        sendSweetAlert(session, "Error", "Main label, Sub label and Feature label must not be identical.", type = "error")
        return()
      }
      
      # Check if the desired columns are present
      desired_cols <- c(Main_label, Sub_label, Feature_label)
      if (!all(desired_cols %in% colnames(data))) {
        sendSweetAlert(session, "Error", "The data does not contain all the necessary columns.", type = "error")
        return()
      }
      
      
      sample_cols <- rownames(seq)
      # Dynamically select columns using the input strings
      df_for_or <- data %>%
        select(!!sym(Main_label), !!sym(Sub_label), !!sym(Feature_label), all_of(sample_cols))
      
      # Remove rows where the Main_label column does not have the Selected_Groups 
      df_for_or <- df_for_or %>% filter(!!sym(Main_label) %in% Selected_Groups)
      
      print(head(df_for_or))
      
      # Remove duplicate rows based on the feature column
      # df_for_or <- df_for_or[!duplicated(df_for_or[[Feature_label]]), ]
      
      # Make the feature names unique 
      df_for_or[[Feature_label]] <- make.unique(as.character(df_for_or[[Feature_label]]))
      
      # Pivot the data from wide to long format using the sample columns
      df_long <- df_for_or %>%
        pivot_longer(
          cols = all_of(sample_cols),
          names_to = "sample",
          values_to = "intensity"
        )
      
      # Make the sequence data frame available for joining by converting rownames to a "sample" column
      seq <- seq %>% rownames_to_column(var = "sample")
      
      # Join the sample information from the sequence to the long-format data
      df_long <- df_long %>%
        left_join(seq, by = "sample") %>%
        select(-batch, -order, -time, -paired, -amount)
      
      # Remove rows with NA intensity and apply log2 transformation
      df_long <- df_long %>% filter(!is.na(intensity))
      # df_long$intensity <- log2(df_long$intensity)
      df_long <- df_long %>% filter(!is.na(intensity), !is.infinite(intensity))
      
      # Filter the long data to only include the two groups of interest
      df_long <- df_long %>%
        filter(group %in% c(group1, group2)) %>%
        mutate(group = factor(group, levels = c(group1, group2)))
      
      # Calculate odds ratios for each feature (using the feature column provided by input)
      odds_ratios <- df_long %>%
        group_by(feature = .data[[Feature_label]]) %>%
        do({
          mod <- glm(group ~ intensity, data = ., family = binomial) # Logistic regression
          broom::tidy(mod) # Extract coefficients
        }) %>%
        filter(term == "intensity") %>%
        mutate(odds_ratio = estimate,
               lower_ci   = (estimate - 1.96 * std.error),
               upper_ci   = (estimate + 1.96 * std.error)) #%>% 
        # mutate(odds_ratio = exp(estimate),
        #        lower_ci   = exp(estimate - 1.96 * std.error),
        #        upper_ci   = exp(estimate + 1.96 * std.error))
      
      # Join the odds ratios back with the original features for label info
      odds_ratios <- odds_ratios %>%
        left_join(df_for_or, by = c("feature" = Feature_label)) %>%
        select(feature, odds_ratio, lower_ci, upper_ci,
               main = !!sym(Main_label),
               sub  = !!sym(Sub_label))
      
      # Determine significance based on the confidence interval
      odds_ratios <- odds_ratios %>%
        mutate(Significance = case_when(
          lower_ci > 1 ~ "Positive",
          upper_ci < 1 ~ "Negative",
          TRUE         ~ "Not significant"
        ))
      
      print(head(odds_ratios))
      message("Max odds ratio:")
      print(max(odds_ratios$odds_ratio))
      message("Min CI:")
      print(min(odds_ratios$lower_ci))
      message("Max CI:")
      print(max(odds_ratios$upper_ci))
      
      # Render the plot
      output$OR_plot <- renderPlot({
        ggplot(odds_ratios, aes(x = odds_ratio, y = sub, color = Significance)) +
          geom_point(position = position_dodge(width = 0.5), size = 4) +
          geom_quasirandom() +
          geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci),
                         height = 0, alpha = 0.5,
                         position = position_dodge(width = 0.5),
                         linetype = 3, linewidth = 2) +
          geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
          scale_x_log10() +
          scale_color_manual(values = c("Negative" = "#1f78b4",
                                        "Not significant" = "gray60",
                                        "Positive" = "#e31a1c")) +
          # Facet by the main label (e.g. super_class or main_class) – labels are now on the left
          facet_grid(main ~ ., scales = "free_y", space = "free_y", switch = "y") +
          theme_bw(base_size = 12) +
          theme(
            strip.placement = "outside",
            strip.text.y.left = element_text(angle = 90),
            legend.position = "right"
          ) +
          labs(
            x = "Odds Ratio (95% CI)",
            y = NULL,
            color = "Direction",
            title = "Odds Ratios by Lipid Class"
          )
      })
      
      # Create a table of significant lipids
      sig_results <- odds_ratios %>% filter(Significance != "Not significant")
      sig_results <- sig_results %>%
        relocate(c(odds_ratio, lower_ci, upper_ci), .after = sub)
      
      output$OR_table <- DT::renderDataTable({
        DT::datatable(sig_results,
                  options = list(
                    scrollX = TRUE,
                    pageLength = 20
                  ))
      })
      
      message(sample(quotes, 1))
    }
  })
  
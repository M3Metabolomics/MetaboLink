  # Update selected data
  observeEvent(input$selectDataset, ignoreInit = TRUE, {
    rv$activeFile <- which(rv$choices %in% input$selectDataset)
    
    output$seq_table <- renderDT(rv$sequence[[rv$activeFile]], extensions = 'Responsive', server = F, 
                                 editable = T, selection = 'none', options = list(pageLength = nrow(rv$sequence[[rv$activeFile]]), 
                                                                                  scrollX = TRUE))
    
    # Make a debugging statement for the active file
    print(paste0("Active file is ", rv$activeFile))
    # Make a debugging statement for the choices
    print(paste0("Choices are ", rv$choices))
    # Make a debugging statement for the selectDataset
    print(paste0("Select Dataset is ", input$selectDataset))
    
    output$seq_table <- renderDT(rv$sequence[[rv$activeFile]],
                                 extensions = 'Responsive',
                                 server = F,
                                 editable = T,
                                 selection = 'none',
                                 options = list(pageLength = nrow(rv$sequence[[rv$activeFile]]),
                                                scrollX = TRUE))
    output$diboxtitle <- renderText(names(rv$data[rv$activeFile]))
    
    output$dttable <- renderDT(rv$data[[rv$activeFile]], rownames = FALSE, options = list(scrollX = TRUE, scrollY = "700px"))
    
    output$dt_drift_panel <- renderDT(rv$data[[rv$activeFile]][rv$sequence[[rv$activeFile]][, 1] %in% "Name"], rownames = FALSE, 
                                      options = list(autoWidth = TRUE, scrollY = "700px", pageLength = 20))
    
    output$dt_boxplot_panel <- renderDT(rv$data[[rv$activeFile]][rv$sequence[[rv$activeFile]][, 1] %in% "Name"],
                                        rownames = FALSE, 
                                        options = list(autoWidth = TRUE,
                                                       scrollY = "700px",
                                                       pageLength = 20))
    
    data <- rv$data[[rv$activeFile]]
    sequence <- rv$sequence[[rv$activeFile]]
    
    output$histogram <- renderPlotly({
      samples <- data[, sequence[ , 'labels'] %in% "Sample"]
      medians <- apply(samples, 2, median, na.rm = TRUE)
      median_data <- data.frame(
        Sample = names(medians),
        Median = medians
      )
      ggplot(median_data, aes(x = Sample, y = Median)) +
        geom_col(fill = "skyblue", color = "black", width = 0.7) +
        labs(x = "Samples", y = "Median") +
        theme_minimal() +
        theme(axis.text.x = element_blank())
    })
    
    output$histogram_qc <- renderUI({
      QCs <- data[, sequence[ , 'labels'] %in% "QC"]
      if(ncol(QCs) > 0) {
        medians <- apply(QCs, 2, median, na.rm = TRUE)
        median_QC <- data.frame(
          QC = names(medians),
          Median = medians
        )
        plotlyOutput("qc_distribution")
        output$qc_distribution <- renderPlotly({
          ggplot(median_QC, aes(x = QC, y = Median)) +
            geom_col(fill = "skyblue", color = "black") +
            labs(x = "Sample", y = "Median") +
            theme_minimal()
        })
      }
      else {
        textOutput("No columns labeled QC.")
      } 
    })
    output$dt_boxplot_panel <- renderDT(rv$data[[rv$activeFile]][rv$sequence[[rv$activeFile]][, 1] %in% "Name"],
                                        rownames = FALSE, 
                                        options = list(autoWidth = TRUE,
                                                       scrollY = "700px",
                                                       pageLength = 20))
    
    
    # Plot of classes in the identifer table
    output$class_plot <- renderPlotly({
      req(input$selected_column_class_plot)  # Ensure a column is selected
      
      # Copy dataset
      data_sub <- data
      colname <- input$selected_column_class_plot
      
      # calculate how many invalid values are in the selected column
      invalid_values <- sum(is.na(data_sub[[colname]]) | data_sub[[colname]] == "" | 
                              data_sub[[colname]] == " " | data_sub[[colname]] == "NA" | 
                              data_sub[[colname]] == "N/A")
      
      # Remove rows where the selected column is NA, empty, or contains unwanted strings.
      data_sub <- data_sub[!is.na(data_sub[[colname]]) &
                             data_sub[[colname]] != "" &
                             data_sub[[colname]] != " " &
                             data_sub[[colname]] != "NA" &
                             data_sub[[colname]] != "N/A", ]
      
      # Check if the column exists (should always be true if choices are set correctly)
      if (!(colname %in% colnames(data_sub))) {
        showNotification(
          paste("No", colname, "column found in the dataset. Make sure a column named", colname, "is present by running 'Gather Identifiers'."),
          type = "message",
          duration = 10
        )
        return(NULL)
      }
      
      # Create the ggplot, sorting and converting the selected column to a factor
      p <- data_sub %>%
        arrange(.data[[colname]]) %>%
        mutate(!!colname := factor(.data[[colname]], levels = unique(.data[[colname]]))) %>%
        ggplot(aes(x = !!sym(colname))) +
        geom_bar(fill = "steelblue", color = "black", width = 0.4) +
        geom_text(stat = "count", aes(label = after_stat(count)), vjust = -15, size = 4.5, color = "red") +
        labs(
          title = paste0("Distribution of ", colname, " - number of features: ", nrow(data), 
                         "\n Invalid values: ", invalid_values),
          x = colname,
          y = "Count"
        ) +
        theme_bw(base_size = 11) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      # Convert the ggplot object to an interactive Plotly plot
      ggplotly(p)
    })
    
    if (sum(rv$sequence[[rv$activeFile]][, 1] %in% "Name") == 1) { #TODO check if this check is needed
      internalStandards <- findInternalStandards(rv$data[[rv$activeFile]][rv$sequence[[rv$activeFile]][, 1] %in% "Name"])
      updateCheckboxGroupInput(session, "isChoose", choices = internalStandards, selected = internalStandards)
      enable("normalizeIS"); enable("removeIS"); enable("saveIS")
      if(length(internalStandards) == 0) {
        disable("normalizeIS"); disable("removeIS"); disable("saveIS")
      } 
    }
  })

  observe({
    req(rv$activeFile,
        rv$data[[rv$activeFile]],
        rv$sequence[[rv$activeFile]])
    
    data <- rv$data[[rv$activeFile]]
    sequence <- rv$sequence[[rv$activeFile]]
    
    output$histogram <- renderPlotly({
      
      message("Inside group histogram")
      # Filter sample columns based on sequence labels
      samples <- data[, sequence[, "labels"] %in% "Sample"]
      
      # Compute medians for each sample
      medians <- apply(samples, 2, median, na.rm = TRUE)
      median_data <- data.frame(
        Sample = names(medians),
        Median = medians,
        stringsAsFactors = FALSE
      )
      
      # Add group information if present
      if ("group" %in% colnames(sequence) && !all(is.na(sequence[, "group"]))) {
        group_data <- data.frame(
          Sample = colnames(data),
          Group = sequence[, "group"],
          stringsAsFactors = FALSE
        )
        median_data <- merge(median_data, group_data, by = "Sample", all.x = TRUE)
        
        # Ensure the samples are ordered within each group
        median_data <- median_data %>%
          arrange(Group, mixedorder(Sample)) %>% # Use mixedorder for natural sorting
          mutate(
            Group = factor(Group, levels = unique(Group)), # Preserve group order
            Sample = factor(Sample, levels = unique(Sample)) # Preserve sample order
          )
        
        # Create a grouped bar chart with distinct group colors
        p <- ggplot(median_data,
                    aes(x = Sample,
                        y = Median,
                        fill = Group)) +
          geom_col(color = "black", width = 0.7) +
          labs(x = "Samples", y = "Median", fill = "Group") +
          facet_wrap(~Group, scales = "free_x", nrow = 1) + 
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            strip.text = element_text(size = 12, face = "bold")
          )
      } else {
        # Ensure samples are sorted naturally (numeric first, then lexicographical)
        median_data <- median_data %>%
          arrange(mixedorder(Sample)) %>% # Use mixedorder for natural sorting
          mutate(Sample = factor(Sample, levels = unique(Sample)))
        
        p <- ggplot(median_data, aes(x = Sample, y = Median)) +
          geom_col(fill = "skyblue", color = "black", width = 0.7) +
          labs(x = "Samples", y = "Median") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
      
      ggplotly(p)
    })
    
    output$histogram_qc <- renderUI({
      QCs <- data[, sequence[, "labels"] %in% "QC"]
      if (ncol(QCs) > 0) {
        medians <- apply(QCs, 2, median, na.rm = TRUE)
        median_QC <- data.frame(
          QC = names(medians),
          Median = medians
        )
        plotlyOutput("qc_distribution")
        output$qc_distribution <- renderPlotly({
          ggplot(median_QC, aes(x = QC, y = Median)) +
            geom_col(fill = "skyblue", color = "black") +
            labs(x = "QC", y = "Median") +
            theme_minimal()
        })
      } else {
        textOutput("No columns labeled QC.")
      }
    })
    
    output$violin_plot <- renderPlotly({
      message("Inside violin plot")
      
      seq_sub <- sequence[sequence$labels == "Sample", ]

      data_sub <- data[, rownames(seq_sub)] 
      
      seq_sub <- seq_sub %>% 
        rownames_to_column("sample")
      
      data_long <- data_sub %>%
        pivot_longer(cols = everything(), names_to = "sample", values_to = "value")
      
      data_long <- left_join(data_long, seq_sub, by = "sample")
      
      data_long <- data_long[!is.na(data_long$value), ]
      
      sample_size <- data_long %>%
        group_by(group) %>%
        summarize(num = n())
      
      data_long <- data_long %>%
        left_join(sample_size, by = "group") %>%
        mutate(myaxis = paste0(group, "\n", "n=", num))
      
      print(head(data_long))
      print(head(sample_size))
      
      # Create the base plot
      violin_layer <- ggplot(data_long, aes(x = myaxis, y = value, fill = group)) +
        geom_violin(width = 1.4) + 
        # scale_fill_viridis(discrete = TRUE) +
        theme_bw() +
        theme(
          legend.position = "none",
          plot.title = element_text(size = 11)
        ) +
        ggtitle("Violin Plot with Boxplot and Sample Size") +
        xlab("")
      
      boxplot_layer <- geom_boxplot(width = 0.1, color = "red", alpha = 0.2)
      
      # Combine them into the final plot
      violin_plot <- violin_layer + boxplot_layer
      
      # For interactivity with plotly:
      ggplotly(violin_plot)
    })
    
  })

  observeEvent(rv$choices, {
    choices <- rv$choices
    num_datasets <- length(choices)
    
    output$downloadSequence <- downloadHandler(
      filename <- function() {
        paste0(names(rv$data[rv$activeFile]), "_seq.csv")
      },
      content = function(file) {
        write.csv(cbind("sample" = rownames(rv$sequence[[rv$activeFile]]), rv$sequence[[rv$activeFile]]), file, row.names = FALSE) # TODO
      }
    )
    
    # Export panel
    output$export_ui <- renderDownloadUI("dwn_general", ".csv")
    output$export_metabo <- renderDownloadUI("dwn_metabo", "_metabo.csv")
    output$export_stats <- renderDownloadUI("dwn_stats", "_results.xlsx")
    output$export_settings <- renderDownloadUI("dwn_settings", ".txt")
    
    output$downloadLipids <- downloadHandler(
      filename = function() {
        paste0("cleaned_lipid_names_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(rv$multipleLipidNamesDf, file, row.names = FALSE)
      }
    )
    
    #    lapply(seq_len(num_datasets), createDownloadHandler("general", ".csv", rv$data[[rv$activeFile]]))
    #    lapply(seq_len(num_datasets), createDownloadHandler("stats", ".xlsx", rv$results[[rv$activeFile]]))
    #    lapply(seq_len(num_datasets), createDownloadHandler("settings", ".txt", rv$info[[rv$activeFile]]))
    #    lapply(seq_len(num_datasets), createDownloadHandler("metabo", ".csv", getMetaboData))
    
    lapply(1:length(rv$choices), function(x) {
      output[[paste0("dwn_stats", x)]] <- downloadHandler(
        filename = function() {
          paste0(names(rv$data[x]), "_results.xlsx")
        },
        content = function(file) {
          write_xlsx(rv$results[[x]], file)
        }
      )
    })
    lapply(1:length(rv$choices), function(x) {
      output[[paste0("dwn_general", x)]] <- downloadHandler(
        filename = function() {
          paste0(names(rv$data[x]), ".csv")
        },
        content = function(file) {
          write.csv(rv$data[[x]], file, row.names = FALSE)
        }
      )
    })
    lapply(1:length(rv$choices), function(x) {
      output[[paste0("dwn_settings", x)]] <- downloadHandler(
        filename = function() {
          paste0(names(rv$data[x]), ".txt")
        },
        content = function(file) {
          write.csv(rv$info[x], file, row.names = FALSE)
        }
      )
    })
    lapply(1:length(rv$choices), function(x) {
      dat <- rv$data[[x]]
      seq <- rv$sequence[[x]]
      seq[seq[, 1] %in% "QC", 4] <- "QC"
      group <- c("", seq[seq[, 1] %in% c("Sample", "QC"), 4])
      outdat <- data.frame(dat[seq[, 1] %in% "Name"], dat[seq[, 1] %in% c("Sample", "QC")])
      outdat <- rbind(group, outdat)
      output[[paste0("dwn_metabo", x)]] <- downloadHandler(
        filename = function() {
          paste0(names(rv$data[x]), "_metabo.csv")
        },
        content = function(file) {
          write.csv(outdat, file, row.names = FALSE)
        }
      )
    })
    
    updateCheckboxGroupInput(session, "export_xml_list", choices = choices, selected = NULL)
    updateCheckboxGroupInput(session, "filesToRemove", choices = names(rv$data), selected = NULL)
    updateSelectInput(session, "drift_select", choices = c("None", choices))
    
    inputs <- c("selectDataset", "mergeFile",
                "selectpca1", "selectpca2",
                "select_data_for_enrichment",
                "select_data_circular_barplot",
                "select_heatmap_data", "select_volcano_data",
                "select_OR_data") # Update select inputs
    
    selected <- ifelse(is.null(rv$activeFile), length(choices), rv$activeFile)
    for(input in inputs) {
      updateSelectInput(session, input, choices = choices, selected = choices[selected])
    }
    
  })

  observeEvent(input$removeFiles, {
    if (is.null(input$filesToRemove)) {
      showNotification("No files selected.", type = "error")
    } else if(length(input$filesToRemove) == length(rv$choices)) {
      showNotification("At least one file must be kept.", type = "error")
    } else {
      keep <- !names(rv$data) %in% input$filesToRemove      
      rv$data <- rv$data[keep]
      rv$sequence <- rv$sequence[keep]
      rv$info <- rv$info[keep]
      rv$results <- rv$results[keep]
      rv$activeFile <- length(rv$data)
      rv$choices <- paste(seq_along(rv$data), ": ", names(rv$data))
      showNotification("Files removed.", type = "message")
    }
  })

  observeEvent(input$reuseSequence, {
    inputSequence <- read.csv(input$inputSequence$datapath, header = 1, stringsAsFactors = FALSE)
    colnames(inputSequence) <- tolower(colnames(inputSequence))
    inputSequence <- checkSequence(inputSequence)
    sequence <- rv$sequence[[rv$activeFile]]
    labeledSequence <- data.frame("sample" = row.names(sequence), sequence)
    inputSequence["sample"] <- lapply(inputSequence["sample"], as.character)
    sequence <- left_join(labeledSequence[, 1:2], inputSequence, by = "sample")
    row.names(sequence) <- sequence[, 1]
    sequence <- sequence[, -1]
    rv$sequence[[rv$activeFile]] <- sequence
  })
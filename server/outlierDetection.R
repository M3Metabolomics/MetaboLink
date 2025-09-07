  # PCA results update 
  # Whenever pca_results are updated, update the selectInput choices for outlier detection
  observe({
    req(rv$pca_results)  # Ensure rv$pca_results is not NULL
    
    input_ids <- c("kmeans_pca", "hierarchical_pca", "dbscan_pca",
                   "hdbscan_pca", "optics_pca", "lof_pca")
    
    choices <- names(rv$pca_results)
    
    # If choices are NULL or empty, provide a default message
    if (is.null(choices) || length(choices) == 0) {
      choices <- list("No PCA results available" = "")
    }
    
    lapply(input_ids, function(input_id) {
      updateSelectInput(session, input_id, choices = choices)
    })
  })
  
  # Function to generate selectInput UI for a given input ID
  create_group_selection_ui <- function(input_id) {
    renderUI({
      if (is.null(rv$activeFile)) {
        showNotification("No data", type = "error")
        return(NULL)
      }
      
      data <- rv$data[[rv$activeFile]]
      sequence <- rv$sequence[[rv$activeFile]]
      seq <- sequence[sequence$labels %in% c("Sample", "QC"), ]
      data <- data[, rownames(seq), drop = FALSE]
      
      if (input[[paste0("select_groups_", input_id)]]) {  # Check if the checkbox is selected
        selectInput(
          inputId = paste0("selected_groups_", input_id),  # Unique inputId per method
          label = paste("Select Groups for", input_id),
          choices = seq$group, 
          selected = seq$group[1],  # Default to first group
          multiple = TRUE, 
          width = "100%"
        )
      }
    })
  }
  
  # Create separate UI outputs for each tab
  output$group_selection_ui_kmeans <- create_group_selection_ui("kmeans")
  output$group_selection_ui_hierarchical <- create_group_selection_ui("hierarchical")
  output$group_selection_ui_dbscan <- create_group_selection_ui("dbscan")
  output$group_selection_ui_hdbscan <- create_group_selection_ui("hdbscan")
  output$group_selection_ui_optics <- create_group_selection_ui("optics")
  output$group_selection_ui_lof <- create_group_selection_ui("lof")
  
  # kmeans
  observeEvent(input$compute_kmeans_eval, {
    # Fetch the PCA data from the reactive values
    pca_data <- rv$pca_results[[input$kmeans_pca]]  # Fetch the selected PCA result (list)
    
    # Check if PCA data is available and extract it
    if (!is.null(pca_data)) {
      pca_df <- pca_data[[1]]  # Extract the pca_df directly from the list
      
      # Get the evaluation method selected by the user
      method <- input$kmeans_eval_method
      
      # Check that the method is valid
      if (!is.null(method)) {
        # Use the selected evaluation method to create the plot
        eval_plot <- switch(method,
                            "silhouette" = fviz_nbclust(pca_df[, c("PC1", "PC2")], kmeans, method = "silhouette", k.max = nrow(pca_df) - 1) + 
                              scale_x_discrete(breaks = seq(1, nrow(pca_df) - 1, by = 5)) +  # Reduce number of x-axis ticks
                              labs(title = "Optimal number of clusters (Silhouette Method)") +
                              theme_bw(),
                            "wss" = fviz_nbclust(pca_df[, c("PC1", "PC2")], kmeans, method = "wss", k.max = nrow(pca_df) - 1) +
                              scale_x_discrete(breaks = seq(1, nrow(pca_df) - 1, by = 5)) +  # Reduce number of x-axis ticks
                              labs(title = "Optimal number of clusters (Elbow Method)") +
                              theme_bw(),
                            "gap_stat" = fviz_nbclust(pca_df[, c("PC1", "PC2")], kmeans, method = "gap_stat", k.max = nrow(pca_df) - 1) +
                              scale_x_discrete(breaks = seq(1, nrow(pca_df) - 1, by = 5)) +  # Reduce number of x-axis ticks
                              labs(title = "Optimal number of clusters (Gap Statistic Method)") +
                              theme_bw())
        
        # Render the evaluation plot in plotly
        output$kmeans_eval_plot <- renderPlotly(ggplotly(eval_plot))
      }
    }
  })
  observeEvent(input$run_kmeans, {
    # Fetch the selected PCA data
    pca_data <- rv$pca_results[[input$kmeans_pca]]  # Fetch the selected PCA result (list)
    
    if (!is.null(pca_data)) {
      pca_df <- pca_data[[1]]   # Extract pca_df from the list
      PC_df <- pca_data[[2]]    # Extract PC_df from the list
      k <- input$num_clusters   # Get the number of clusters (k)
      percentile_threshold <- input$percentile_threshold  # Get the percentile threshold
      
      # in pca_df remove rows with QC in column group 
      pca_df <- pca_df[pca_df$group != "QC", ]
      
      print(pca_df)
      print(PC_df)
      pca_df$group <- as.character(pca_df$group)
      
      if (input$select_groups_kmeans) {
        message(paste0("Enable grouping for groups: ", input$select_groups_kmeans))
        if (is.null(input$selected_groups_kmeans) || length(input$selected_groups_kmeans) < 1) {
          showNotification("Select at least one groups for the heatmap.", type = "error")
          return()  # Stop execution
        }
        # Filter seq_subset and data_subset by selected groups
        selected_groups_kmeans <- input$selected_groups_kmeans
        message(paste0("Selected groups: ", paste(selected_groups_kmeans, collapse = ", ")))
        pca_df <- pca_df[pca_df$group %in% selected_groups_kmeans, ]
      }
      
      pca_df$group <- as.factor(pca_df$group)
      
      print(head(pca_df))
      print(head(PC_df))

      # Check if pca_df is a valid data frame
      if (!is.null(pca_df) && is.data.frame(pca_df)) {
        # Call the kmeans_clustering function to get the results
        kmeans_results <- kmeans_clustering(pca_df, k, percentile_threshold, PC_df)
        
        # Generate a unique name for the K-means result from the selected PCA result name
        pca_result_name <- input$kmeans_pca   # This gives the selected name from the dropdown
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")  # Timestamp for uniqueness (date + time)
        kmeans_name <- paste0(pca_result_name, "_Kmeans_", timestamp)
        
        # Save the K-means result (including outlier information) into the reactive value
        rv$outlier_df[[kmeans_name]] <- list(kmeans_df = kmeans_results)
        
        # Render the K-means plot in the UI
        output$kmeans_plot <- renderPlotly({
          kmeans_results$kmeans_plotly  # Pass the interactive plotly object to renderPlotly
        })
        
        # Render the K-means outlier table in the UI
        output$kmeans_outliers <- renderDT({
          datatable(kmeans_results$kmeans_df, options = list(pageLength = 20, autoWidth = TRUE))
        })
        
      } else {
        cat("Selected PCA result is not a valid data frame.\n")
      }
    } else {
      cat("No PCA results selected or available.\n")
    }
  })
  # Hierarchical Clustering 
  observeEvent(input$run_hierarchical, {
    # Fetch the selected PCA data from reactive values
    pca_data <- rv$pca_results[[input$hierarchical_pca]]  # Fetch the selected PCA result (list)
    sequence <- rv$sequence[[rv$activeFile]]  # Fetch the sequence data from the active file
    
    seq_subset <- sequence[sequence[, "labels"] %in% c("Sample", "QC"), ]  # Get the sequence for the samples and QC
    
    # Debugging to show the sequence data
    if (!"labels" %in% colnames(sequence)) {
      cat("Error: 'labels' column not found in sequence data.\n")
      return()
    }
    str(seq_subset)
    
    if (!is.null(pca_data)) {
      pca_df <- pca_data[[1]]   # Extract pca_df from the list
      PC_df <- pca_data[[2]]    # Extract PC_df from the list
      
      # in pca_df remove rows with QC in column group 
      pca_df <- pca_df[pca_df$group != "QC", ]
      
      # Retrieve parameters from the UI
      method <- input$clustering_method   # Clustering method selected
      k <- input$num_clusters_hierarchical  # Number of clusters
      threshold <- input$threshold  # Dendrogram threshold
      
      pca_df$group <- as.character(pca_df$group)
      
      if (input$select_groups_hierarchical) {
        message(paste0("Enable grouping for groups: ", input$select_groups_hierarchical))
        if (is.null(input$selected_groups_hierarchical) || length(input$selected_groups_hierarchical) < 1) {
          showNotification("Select at least one groups for the hierarchical.", type = "error")
          return()  # Stop execution
        }
        # Filter seq_subset and data_subset by selected groups
        selected_groups_hierarchical <- input$selected_groups_hierarchical
        message(paste0("Selected groups: ", paste(selected_groups_hierarchical, collapse = ", ")))
        pca_df <- pca_df[pca_df$group %in% selected_groups_hierarchical, ]
      }
      
      pca_df$group <- as.factor(pca_df$group)
      
      # Debugging to show the parameters 
      # cat("Clustering method: ", method, "\n")
      # cat("Number of clusters: ", k, "\n")
      # cat("Dendrogram threshold: ", threshold, "\n")
      # print(head(seq_subset))
      
      # Perform the hierarchical clustering
      hierarchical_results <- perform_hierarchical_clustering(pca_df, seq_subset, method, k, threshold)
      
      # Render the hierarchical clustering plot
      output$hclust_plot <- renderPlotly({
        hierarchical_results$hclust_plot  # Pass the plotly object to renderPlotly
      })
      
      # conf_matrix_plot <- perform_confusion_matrix(pca_df, seq_subset, method, k)
      
      # Render the confusion matrix plot
      # output$conf_matrix_plot <- renderPlotly({
      #   conf_matrix_plot  # Pass the plotly object to renderPlotly
      # })
      
      # Render the dendrogram plot
      output$dendrogram_plot <- renderPlotly({
        hierarchical_results$dendrogram_plot  # Pass the plotly dendrogram to renderPlotly
      })
      
      # Render the outliers table
      output$hierarchical_outliers <- renderDT({
        datatable(hierarchical_results$hierarchical_outliers, options = list(pageLength = 20, autoWidth = TRUE))
      })
      
    } else {
      cat("No PCA results selected or available.\n")
    }
  })

    # Run DBSCAN
  # DBSCAN Clustering
  observeEvent(input$compute_knn, {
    req(rv$pca_results[[input$dbscan_pca]])  # Ensure PCA data is loaded
    pca_data <- rv$pca_results[[input$dbscan_pca]]  # Fetch the selected PCA result (list)
    pca_df <- pca_data[[1]]  # Extract pca_df from the list
    k <- input$knn
    
    # in pca_df remove rows with QC in column group 
    pca_df <- pca_df[pca_df$group != "QC", ]
    
    
    # Debug print
    # print(paste("Computing kNN distance plot with k =", k))
    
    output$knn_plot <- renderPlotly({
      perform_kNN_dist_plot(pca_df, k)
    })
  })
  observeEvent(input$run_dbscan, {
    req(rv$pca_results[[input$dbscan_pca]])  # Ensure PCA data is loaded
    pca_data <- rv$pca_results[[input$dbscan_pca]]  # Fetch the selected PCA result (list)
    pca_df <- pca_data[[1]]  # Extract pca_df from the list
    
    # in pca_df remove rows with QC in column group 
    pca_df <- pca_df[pca_df$group != "QC", ]
    
    pca_df$group <- as.character(pca_df$group)
    
    if (input$select_groups_hdbscan) {
      message(paste0("Enable grouping for groups: ", input$select_groups_hdbscan))
      if (is.null(input$selected_groups_hdbscan) || length(input$selected_groups_hdbscan) < 1) {
        showNotification("Select at least one groups for the hdbscan.", type = "error")
        return()  # Stop execution
      }
      # Filter seq_subset and data_subset by selected groups
      selected_groups_hdbscan <- input$selected_groups_hdbscan
      message(paste0("Selected groups: ", paste(selected_groups_hdbscan, collapse = ", ")))
      pca_df <- pca_df[pca_df$group %in% selected_groups_hdbscan, ]
    }
    
    pca_df$group <- as.factor(pca_df$group)

    eps <- input$eps
    min_pts <- input$min_pts_dbscan
    
    # Debugging print
    # print(paste("Running DBSCAN with eps =", eps, "and minPts =", min_pts))
    
    dbscan_res <- perform_dbscan_clustering(pca_df, eps, min_pts)
    
    dbscan_res$dbscan_outliers <- dbscan_res$dbscan_outliers %>%
      rename(Category = Outlier)
    
    output$dbscan_plot <- renderPlotly({
      dbscan_res$dbscan_plot
    })
    # Debugging 
    # print(head(dbscan_res$dbscan_outliers))
    # print(names(dbscan_res$dbscan_outliers))
    
    output$dbscan_outliers <- renderDT({
      datatable(dbscan_res$dbscan_outliers %>%
                  select(Sample, PC1, PC2, Cluster, Category), options = list(pageLength = 10))
    })
  })


    # HDBSCAN Clustering
  observeEvent(input$run_hdbscan, {
    # Ensure PCA data is loaded
    req(rv$pca_results[[input$hdbscan_pca]])
    
    # Fetch the selected PCA result
    pca_data <- rv$pca_results[[input$hdbscan_pca]]
    pca_df <- pca_data[[1]]  # Extract pca_df from the list
    
    # in pca_df remove rows with QC in column group 
    pca_df <- pca_df[pca_df$group != "QC", ]
    
    min_pts <- input$min_pts_hdbscan  # Minimum number of points for HDBSCAN
    threshold <- input$threshold_hdbscan  # Outlier threshold
    
    pca_df$group <- as.character(pca_df$group)
    
    if (input$select_groups_hdbscan) {
      message(paste0("Enable grouping for groups: ", input$select_groups_hdbscan))
      if (is.null(input$selected_groups_hdbscan) || length(input$selected_groups_hdbscan) < 1) {
        showNotification("Select at least one groups for the hdbscan.", type = "error")
        return()  # Stop execution
      }
      # Filter seq_subset and data_subset by selected groups
      selected_groups_hdbscan <- input$selected_groups_hdbscan
      message(paste0("Selected groups: ", paste(selected_groups_hdbscan, collapse = ", ")))
      pca_df <- pca_df[pca_df$group %in% selected_groups_hdbscan, ]
    }
    
    pca_df$group <- as.factor(pca_df$group)
    
    # Debug print
    # print(paste("Running HDBSCAN with minPts =", min_pts))
    
    # Perform HDBSCAN clustering
    hdbscan_res <- perform_hdbscan_clustering(pca_df, min_pts)
    
    # Process outliers from the HDBSCAN results
    hdbscan_outliers <- hdbscan_res$hdbscan_outliers
    
    # Add categorization based on OutlierScore and threshold
    if (nrow(hdbscan_outliers) > 0 && "OutlierScore" %in% names(hdbscan_outliers)) {
      hdbscan_outliers <- hdbscan_outliers %>%
        mutate(Category = ifelse(OutlierScore > threshold, "Outlier", "Inlier"))
    } else {
      showNotification("HDBSCAN results are empty or do not contain 'OutlierScore' column.", type = "error")
    }
    
    # Render HDBSCAN plot
    output$hdbscan_plot <- renderPlotly({
      hdbscan_res$hdbscan_plot
    })
    
    # Render outlier table
    output$hdbscan_outliers <- renderDT({
      datatable(hdbscan_outliers %>%
                  select(Sample, PC1, PC2, Cluster, Category), options = list(pageLength = 10))
    })
  })
  # OPTICS Clustering
  observeEvent(input$run_optics, {
    req(rv$pca_results[[input$optics_pca]])
    pca_data <- rv$pca_results[[input$optics_pca]]
    pca_df <- pca_data[[1]]  # Extract pca_df from the list
    
    # in pca_df remove rows with QC in column group 
    pca_df <- pca_df[pca_df$group != "QC", ]
    
    # Debugging print
    min_pts <- input$min_pts_optics
    eps <- if (is.na(input$eps_optics)) NULL else input$eps_optics
    eps_cl <- input$eps_cl_optics
    
    pca_df$group <- as.character(pca_df$group)
    
    if (input$select_groups_optics) {
      message(paste0("Enable grouping for groups: ", input$select_groups_optics))
      if (is.null(input$selected_groups_optics) || length(input$selected_groups_optics) < 1) {
        showNotification("Select at least one groups for the optics.", type = "error")
        return()  # Stop execution
      }
      # Filter seq_subset and data_subset by selected groups
      selected_groups_optics <- input$selected_groups_optics
      message(paste0("Selected groups: ", paste(selected_groups_optics, collapse = ", ")))
      pca_df <- pca_df[pca_df$group %in% selected_groups_optics, ]
    }
    
    pca_df$group <- as.factor(pca_df$group)
    
    # print(paste("Running OPTICS with minPts =", min_pts, ", eps =", eps, ", and eps_cl =", eps_cl))
    
    optics_res <- perform_optics_analysis(pca_df, eps, min_pts, eps_cl)
    
    optics_outliers <- optics$optics_outliers
    
    output$optics_reachability_plot <- renderPlot({
      optics_res$reachability_plot()
    })
    
    output$reachability_plot_threshold <- renderPlot({
      optics_res$reachability_plot_threshold()
    })
    
    output$cluster_plot <- renderPlot({
      optics_res$cluster_plot()
    })
    
    output$optics_outliers <- renderTable({
      optics_res$optics_outliers
    })
  })
  # LOF Clustering
  observeEvent(input$run_lof, {
    req(rv$pca_results[[input$optics_pca]])
    pca_data <- rv$pca_results[[input$optics_pca]]
    pca_df <- pca_data[[1]]  # Extract pca_df from the list
    
    # in pca_df remove rows with QC in column group 
    pca_df <- pca_df[pca_df$group != "QC", ]
    
    threshold <- input$lof_threshold
    min_pts <- input$lof_k
    
    pca_df$group <- as.character(pca_df$group)
    
    if (input$select_groups_lof) {
      message(paste0("Enable grouping for groups: ", input$select_groups_lof))
      if (is.null(input$selected_groups_lof) || length(input$selected_groups_lof) < 1) {
        showNotification("Select at least one groups for the lof. ", type = "error")
        return()  # Stop execution
      }
      # Filter seq_subset and data_subset by selected groups
      selected_groups_lof <- input$selected_groups_lof
      message(paste0("Selected groups: ", paste(selected_groups_lof, collapse = ", ")))
      pca_df <- pca_df[pca_df$group %in% selected_groups_lof, ]
    }
    
    pca_df$group <- as.factor(pca_df$group)
    
    lof_res <- calculate_and_plot_lof(pca_df, threshold = threshold, minPts = min_pts)
    
    lof_plot <- lof_res$lof_plotly
    lof_od_plot <- lof_res$lof_od_plotly
    lof_outliers <- lof_res$lof_outliers
    
    output$lof_plot <- renderPlotly({
      lof_res$lof_plotly
    })
    
    output$lof_od_plot <- renderPlotly({
      lof_res$lof_od_plotly
    })
    
    output$lof_outliers <- renderTable({
      lof_res$lof_outliers
    })
  })

  
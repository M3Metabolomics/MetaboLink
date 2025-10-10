mfa <- reactiveValues(
  selected = NULL,
  threshold_data = NULL,
  tracer_data = NULL,
  tracer_sequence = NULL,
  normalized_a0 = NULL,
  normalized_sum = NULL,
  metabolites = NULL,
  isotopologues = NULL,
  samples = NULL,
  groups = NULL,
  group_time = NULL
)

plot_settings <- reactiveValues(
  intensity_threshold = 500,
  metabolite = NULL,
  sample = NULL,
  metabolite_group = NULL,
  group = NULL,
  group_iso = NULL,
  plot_type = "barplot",
  metabolite_time_table = NULL,
  group_time = NULL,
  plot_type_time = "barplot"

)

isotopologue_settings <- reactiveValues(
  metabolite_iso = NULL,
  group_iso = NULL,
  isotopologues = NULL,
  show_top5 = FALSE,
  group_time = NULL,
  data_type = "raw",
  plot_type_isotopologues = "barplot"
)


observeEvent(input$is_tracer_data, {
    #TODO: fill UI options
    if(input$is_tracer_data) {
      mfa$tracer_data <- rv$data[[rv$activeFile]]
      mfa$metabolites <- unique(mfa$tracer_data$Analyte)
      mfa$isotopologues <- grep("^A\\+", colnames(mfa$tracer_data), value = TRUE)
      mfa$samples <- unique(mfa$tracer_data$Analysis)

      updateSelectInput(session, "metabolite", choices = mfa$metabolites)
      updateSelectInput(session, "metabolite_iso", choices = mfa$metabolites)
      updateSelectInput(session, "metabolite_group", choices = mfa$metabolites)
      updateSelectInput(session, "metabolite_time_table", choices = mfa$metabolites)
      updateSelectInput(session, "sample", choices = mfa$samples)

      # Compute normalized datasets used by plots
      # normalized to A.0 (assuming first isotopologue column is A+0 if present)
      iso_cols <- intersect(mfa$isotopologues, colnames(mfa$tracer_data))
      if (length(iso_cols) > 0) {
        # create matrices of only isotopologue columns
        iso_mat <- as.matrix(mfa$tracer_data[, iso_cols, drop = FALSE])
        # normalized to the +0 isotopologue if present (match any label ending with '+0', e.g. M+0, A+0)
        a0_matches <- grep("\\+0$", iso_cols, value = TRUE)
        a0_col <- if (length(a0_matches) > 0) a0_matches[1] else NA
        if (!is.na(a0_col) && a0_col %in% colnames(iso_mat)) {
          # divide each isotopologue column by the +0 column, preserving matrix structure
          mfa$normalized_a0 <- iso_mat / (iso_mat[, a0_col, drop = TRUE])
        } else {
          mfa$normalized_a0 <- NULL
        }
        # normalized to row sums
        row_sums <- rowSums(iso_mat, na.rm = TRUE)
        # avoid division by zero
        row_sums[row_sums == 0] <- NA
        mfa$normalized_sum <- iso_mat / row_sums
      } else {
        mfa$normalized_a0 <- NULL
        mfa$normalized_sum <- NULL
      }

      
    }
})


observeEvent(input$inputTracerSequence, {
  # read uploaded sequence file safely
  tracer_sequence <- tryCatch(
    read.csv(input$inputTracerSequence$datapath, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (is.null(tracer_sequence)) {
    showNotification("Failed to read tracer sequence file. Please upload a valid CSV/TXT.", type = "error")
    return()
  }

  # normalize column names (case-insensitive) and require sample/group/time
  colmap <- tolower(colnames(tracer_sequence))
  required <- c("sample", "group", "time")
  missing_cols <- setdiff(required, intersect(required, colmap))
  if (length(missing_cols) > 0) {
    showNotification(paste0("Tracer sequence is missing required columns: ", paste(missing_cols, collapse = ", ")), type = "error")
    return()
  }

  # rename columns to a consistent casing so server code can access them reliably
  names(tracer_sequence)[which(colmap == "sample")] <- "sample"
  names(tracer_sequence)[which(colmap == "group")] <- "group"
  names(tracer_sequence)[which(colmap == "time")] <- "time"

  mfa$tracer_sequence <- tracer_sequence

  # compute groups/timepoints safely
  mfa$groups <- unique(tracer_sequence$group)
  mfa$time_points <- sort(unique(tracer_sequence$time))
  mfa$group_time <- unique(paste(tracer_sequence$group, tracer_sequence$time, sep = "_"))

  updateSelectInput(session, "group", choices = mfa$groups)
  updateSelectInput(session, "group_iso", choices = mfa$groups)
  updateSelectInput(session, "time_point", choices = mfa$time_points)
  updatePickerInput(session, "group_time", choices =  mfa$group_time)

  # isotopologue choices should be restricted to those present in tracer data
  if (!is.null(mfa$tracer_data) && !is.null(mfa$isotopologues)) {
    iso_choices <- intersect(mfa$isotopologues, colnames(mfa$tracer_data))
  } else {
    iso_choices <- mfa$isotopologues
  }
  updatePickerInput(session, "isotopologues", choices = iso_choices)

  output$tracer_sequence <- renderDT({
    tracer_sequence
  })

})

# serve example tracer sequence CSV from example_files
output$download_tracer_example <- downloadHandler(
  filename = function() { "tracer_sequence_example.csv" },
  content = function(file) {
    file.copy(file.path("example_files", "tracer_sequence_example.csv"), file)
  }
)


### OVERVIEW PANEL ###

observeEvent(input$update_threshold, {
  req(input$intensity_threshold, mfa$tracer_data)
  threshold <- as.numeric(input$intensity_threshold)

  print(paste("Updating intensity threshold to:", threshold))

  if (!is.null(mfa$tracer_data)) {

    filtered_data <- mfa$tracer_data[rowSums(mfa$tracer_data[, mfa$isotopologues], na.rm = TRUE) >= threshold, ]

    output$tracer_table <- renderDT({
      filtered_data
    })

    mfa$threshold_data <- filtered_data
  }
})


### Isotopologue Profiles ###

output$tracer_plot_ref <- renderPlotly({
  req(mfa$tracer_data, input$metabolite, input$sample)

  plot_settings$metabolite <- input$metabolite
  plot_settings$sample <- input$sample

  plotIsotopologueDistA0(mfa$tracer_data, plot_settings)

})


output$tracer_plot_rowsum <- renderPlot({
  req(mfa$normalized_sum)
  req(input$metabolite)
  req(input$sample)

  plot_settings$metabolite <- input$metabolite
  plot_settings$sample <- input$sample

  plotIsotopologueDist(mfa$tracer_data, plot_settings)
})

#output$tracer_table <- renderDT({
#  req(mfa$tracer_data)
#  req(input$metabolite)
#  req(input$sample)
#
#  plot_settings$metabolite <- input$metabolite
#  plot_settings$sample <- input$sample
#
#  metabolite_data <- mfa$tracer_data[mfa$tracer_data$Analyte == #plot_settings$metabolite & mfa$tracer_data$Analysis == plot_settings$sample, ]
#  metabolite_data
#})


### Fractional Contribution ###

output$fc_plot <- renderPlot({
  req(mfa$tracer_data)
  req(mfa$tracer_sequence)
  req(input$metabolite_group)
  req(input$group)
  req(input$plot_type)

  plot_settings$metabolite_group <- input$metabolite_group
  plot_settings$group <- input$group
  plot_settings$plot_type <- input$plot_type

  # Select samples in group
  group_samples <- mfa$tracer_sequence[mfa$tracer_sequence[,'group'] %in% plot_settings$group, 'sample']
  
  plot_settings$sample <- group_samples

  plotFractionalContribution(mfa$tracer_data, mfa$tracer_sequence, plot_settings)

})

output$fc_table <- renderDT({
  req(mfa$tracer_data)
  req(mfa$tracer_sequence)
  req(input$metabolite_group)
  req(input$group)
  req(input$plot_type)

  plot_settings$metabolite_group <- input$metabolite_group
  plot_settings$group <- input$group
  plot_settings$plot_type <- input$plot_type

  # Select samples in group
  group_samples <- mfa$tracer_sequence[mfa$tracer_sequence[,'group'] %in% plot_settings$group, 'sample']
  
  plot_settings$sample <- group_samples

  selectFCtable(mfa$tracer_data, mfa$tracer_sequence, plot_settings)

})


### Isotopologue Timecourse ###

observeEvent(input$update_iso_plot, {
  req(mfa$tracer_data)
  req(mfa$tracer_sequence)
  req(input$metabolite_iso)
  req(input$group_iso)

  isotopologue_settings$metabolite_iso <- input$metabolite_iso
  isotopologue_settings$group_iso <- input$group_iso
  isotopologue_settings$isotopologues <- input$isotopologues
  isotopologue_settings$show_top5 <- input$show_top5
  isotopologue_settings$group_time <- input$group_time
  isotopologue_settings$data_type <- input$data_type
  isotopologue_settings$plot_type_isotopologues <- input$plot_type_isotopologues

  top5_isotopologues <- NULL
  if(isotopologue_settings$show_top5) {
    # Get top 5 isotopologues by mean abundance (only for isotopologues present in data)
    metabolite_data <- mfa$tracer_data[mfa$tracer_data$Analyte == isotopologue_settings$metabolite_iso, ]
    available_iso <- intersect(mfa$isotopologues, colnames(metabolite_data))
    if (length(available_iso) == 0) {
      showNotification("No isotopologue columns found for selected metabolite.", type = "warning")
    } else {
      mean_abundances <- colMeans(metabolite_data[, available_iso, drop = FALSE], na.rm = TRUE)
      top5_isotopologues <- names(sort(mean_abundances, decreasing = TRUE))[seq_len(min(5, length(mean_abundances)))]
      updatePickerInput(session, "isotopologues", selected = top5_isotopologues)
    }
  }

  # Select relevant data: metabolite, sample
  selected_samples <- mfa$tracer_sequence[mfa$tracer_sequence[,'group'] %in% isotopologue_settings$group_iso, c('sample', 'time')]

  # pick source data depending on data_type
  source_data <- NULL
  if (!is.null(isotopologue_settings$data_type) && isotopologue_settings$data_type == 'normalizedRowSums') {
    if (!is.null(mfa$normalized_sum)) {
      # rebuild a data.frame with isotopologue cols, keeping Analyte and Analysis
      iso_cols <- intersect(mfa$isotopologues, colnames(mfa$tracer_data))
      if (length(iso_cols) > 0) {
        # numeric normalized matrix corresponds row-wise to tracer_data rows
        norm_df <- as.data.frame(mfa$normalized_sum)
        colnames(norm_df) <- iso_cols
        source_data <- cbind(mfa$tracer_data[, c('Analyte', 'Analysis'), drop=FALSE], norm_df)
      }
    } else {
      showNotification('Normalized row-sum data not available. Falling back to raw data.', type = 'warning')
      source_data <- mfa$tracer_data
    }
  } else {
    source_data <- mfa$tracer_data
  }

  selected <- source_data[source_data$Analyte == isotopologue_settings$metabolite_iso & source_data$Analysis %in% selected_samples$sample, ]

  # Keep only selected isotopologues (if top5_isotopologues present, else fall back to isotopologues from input)
  chosen_iso <- NULL
  if (!is.null(top5_isotopologues) && length(top5_isotopologues) > 0) {
    chosen_iso <- intersect(top5_isotopologues, colnames(selected))
  } else if (!is.null(isotopologue_settings$isotopologues) && length(isotopologue_settings$isotopologues) > 0) {
    chosen_iso <- intersect(isotopologue_settings$isotopologues, colnames(selected))
  } else {
    chosen_iso <- intersect(mfa$isotopologues, colnames(selected))
  }

  if (length(chosen_iso) == 0) {
    selected <- selected[, c('Analyte', 'Analysis')]
    showNotification('No isotopologue columns selected or available for plotting.', type = 'warning')
  } else {
    selected <- selected[, c('Analyte', 'Analysis', chosen_iso), drop = FALSE]
  }

  output$iso_table <- renderDT({
    selected
  })

  output$iso_plot <- renderPlot({
    plotStackedIsotopologues(selected, isotopologue_settings)
  })

})

output$isotopologue_plot <- renderPlot({
  req(mfa$tracer_data)
  req(mfa$tracer_sequence)
  req(input$metabolite_iso)
  req(input$group_iso)
  req(input$isotopologues)  
  req(input$show_top5)
  req(input$group_time)
  req(input$data_type)
  req(input$plot_type_isotopologues)

  isotopologue_settings$metabolite_iso <- input$metabolite_iso
  isotopologue_settings$group_iso <- input$group_iso
  isotopologue_settings$isotopologues <- input$isotopologues
  isotopologue_settings$show_top5 <- input$show_top5
  isotopologue_settings$group_time <- input$group_time
  isotopologue_settings$data_type <- input$data_type
  isotopologue_settings$plot_type_isotopologues <- input$plot_type_isotopologues

  plotIsotopologueTimeCourse(mfa$tracer_data, mfa$tracer_sequence, isotopologue_settings)

})


### Group x Time ###


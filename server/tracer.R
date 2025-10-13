
#TODO: make all of these reactive like in MFA3
mfa <- reactiveValues(
  selected = NULL,
  threshold_data = NULL,
  tracer_data = NULL,
  tracer_sequence = NULL,
  normalized_sum = NULL,
  ggplot_data = NULL,
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

observeEvent(input$is_tracer_data, {
    #TODO: validate tracer data format
    if(input$is_tracer_data) {
      mfa$tracer_data <- rv$data[[rv$activeFile]]
      mfa$metabolites <- unique(mfa$tracer_data$Analyte)
      mfa$isotopologues <- grep("^A\\+", colnames(mfa$tracer_data), value = TRUE) #TODO: read also M+X, other patterns?
      mfa$samples <- unique(mfa$tracer_data$Analysis)


      #TODO: function for this
      updateSelectInput(session, "metabolite", choices = mfa$metabolites)
      updateSelectInput(session, "metabolite_iso", choices = mfa$metabolites)
      updateSelectInput(session, "metabolite_group", choices = mfa$metabolites)
      updateSelectInput(session, "metabolite_time_table", choices = mfa$metabolites)
      updateSelectInput(session, "sample", choices = mfa$samples)
      updatePickerInput(session, "isotopologues_iso", choices = mfa$isotopologues)


      #TODO: deal with threshold since it has to be applied before normalization


      # Normalization
      numerical_cols <- mfa$tracer_data[, mfa$isotopologues, drop = FALSE]
      normalized_sum <- numerical_cols / rowSums(numerical_cols, na.rm = TRUE)
      mfa$normalized_sum <- cbind(mfa$tracer_data[, c("Analyte", "Analysis")], normalized_sum)      

    }
})


observeEvent(input$inputTracerSequence, {
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

  #TODO: check that samples in tracer_sequence are present in tracer_data

  mfa$tracer_sequence <- tracer_sequence
  mfa$groups <- unique(tracer_sequence$group)
  mfa$time_points <- sort(unique(tracer_sequence$time))
  mfa$group_time <- unique(paste(tracer_sequence$group, tracer_sequence$time, sep = "_"))
  mfa$tracer_sequence$group_time <- paste(tracer_sequence$group, tracer_sequence$time, sep = "_")

  updateSelectInput(session, "group", choices = mfa$groups)
  updateSelectInput(session, "group_iso", choices = mfa$groups)
  updateSelectInput(session, "time_point", choices = mfa$time_points)
  updatePickerInput(session, "group_time", choices =  mfa$group_time)

  output$tracer_sequence <- renderDT({
    mfa$tracer_sequence
  })

  mfa$ggplot_data <- format4ggplot(mfa$normalized_sum, mfa$tracer_sequence, mfa$isotopologues)

  output$ggplotdata <- renderDT({
    mfa$ggplot_data
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

#output$tracer_plot_ref <- renderPlotly({
#  req(mfa$tracer_data, input$metabolite, input$sample)
#
#  plot_settings$metabolite <- input$metabolite
#  plot_settings$sample <- input$sample
#
#  plotIsotopologueDistA0(mfa$tracer_data, plot_settings)
#
#})

output$tracer_plot_rowsum <- renderPlotly({
  req(mfa$ggplot_data, input$metabolite, input$sample)

  plot_settings$metabolite <- input$metabolite
  plot_settings$sample <- input$sample

  plotIsotopologueDist(mfa$ggplot_data, plot_settings)
})


### Fractional Contribution ###

output$fc_plot <- renderPlotly({
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
  req(mfa$tracer_data, mfa$tracer_sequence, input$metabolite_iso, input$group_iso)

  top5_isotopologues <- NULL
  if(input$show_top5) { #TODO: ask Jesper top5 overall or for this metabolite and group?
    # Get top 5 isotopologues by mean abundance (only for isotopologues present in data)
    metabolite_data <- mfa$normalized_sum[mfa$normalized_sum$Analyte == input$metabolite_iso, ]
    available_iso <- intersect(mfa$isotopologues, colnames(metabolite_data))
    if (length(available_iso) == 0) {
      showNotification("No isotopologue columns found for selected metabolite.", type = "warning")
    } else {
      mean_abundances <- colMeans(metabolite_data[, available_iso], na.rm = TRUE)
      top5_isotopologues <- names(sort(mean_abundances, decreasing = TRUE))[seq_len(min(5, length(mean_abundances)))]
      updatePickerInput(session, "isotopologues", selected = top5_isotopologues)
    }
  }

  #TODO update available isotopologues based on selected metabolite
  #TODO normalize based on selected isotopologues (does it make sense? the other ones are still present in the data) - cannot use ggplot data for this?

  data <- mfa$ggplot_data[mfa$ggplot_data$Analyte == input$metabolite_iso, ]
  data <- data[data$Group %in% input$group_iso, ]
  data <- data[data$Isotopologue %in% if (!is.null(input$isotopologues) && length(input$isotopologues) > 0) {
    input$isotopologues
  } else if (!is.null(top5_isotopologues) && length(top5_isotopologues) > 0) {
    top5_isotopologues
  } else {
    mfa$isotopologues
  }, ]

  
  # Group by group_time 
  original_samples <- data %>%
    filter(!is.na(groupTime)) %>%
    filter(groupTime %in% input$group_time) %>%
    pull(Analysis) %>%
    unique()


  data <- data %>%
    select(Analyte, Isotopologue, Abundance, Analysis, groupTime) %>%
    filter(!is.na(Abundance))

  data <- data %>%
    filter(!is.na(groupTime)) %>%
    filter(groupTime %in% input$group_time)

  # Replicates per time point and normalize
  # 1) normalize within each Analysis for given Analyte so isotopologues sum to 1
  # 2) summarize replicates per analyte + group_time + isotopologue 
  summarized <- data %>%
    group_by(Analyte, Analysis) %>%
    mutate( Abundance = Abundance / sum(Abundance, na.rm = TRUE) ) %>%
    ungroup() %>%
    group_by(groupTime, Isotopologue) %>%
    summarise(
      n_replicates = n_distinct(Analysis),
      Abundance = sum(Abundance, na.rm = TRUE) / n_replicates,
      .groups = 'drop'
    ) %>%
    select(-n_replicates)


  output$iso_table2 <- renderDT({
    summarized
  })


  ## Final samples after filtering
  #final_samples <- summarized %>% pull(Analysis) %>% unique()
  #print(final_samples)
  #excluded_samples <- setdiff(original_samples, final_samples)
  #if(length(excluded_samples) > 0) {
  #  output$exclusion_warning <- renderUI({
  #    div(
  #      style = "color: red;",
  #      paste("Warning: The following samples were excluded due to missing group/time information:", paste(excluded_samples, collapse = ", "))
  #    )
  #  })
  #} else {
  #  output$exclusion_warning <- renderUI({ NULL })
  #}


  output$iso_table <- renderDT({
    data
  })

  settings <- list(
    metabolite_iso = NULL,
    group_time = NULL
  )
  settings$metabolite_iso <- input$metabolite_iso
  settings$group_time <- input$group_time

  output$iso_plot <- renderPlotly({
    plotStackedIsotopologues(summarized, settings)
  })

})


### Group x Time ###



#TODO: make all of these reactive like in MFA3
mfa <- reactiveValues(
  raw = NULL,
  threshold = NULL,
  sequence = NULL,
  normalized_sum = NULL,
  normalized_zero = NULL,
  normalized_long_format = NULL,
  normalized_long_format_ref = NULL
)

metadata <- reactiveValues(
  metabolites = NULL,
  isotopologues = NULL,
  samples = NULL,
  groups = NULL,
  group_time = NULL,
  time_points = NULL
)

# Threshold data - separate for easier access and to keep raw intact
threshold <- reactiveValues(
  data = NULL,
  normalized_sum = NULL,
  normalized_zero = NULL,
  long_format_sum = NULL,
  long_format_zero = NULL
)

threshold_metadata <- reactiveValues(
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
  plot_type_time = "barplot",
  data_type = "normalized_sum",
)

validate_tracer_data <- function(data) {
  required_cols <- c("Analyte", "Analysis")
  missing_required <- setdiff(required_cols, colnames(data))

  isotopologue_pattern <- "^(?:A|M)\\+\\d+$"  # Pattern for isotopologues like A+0, A+1, A+2, etc.
  isotopologues <- grep(isotopologue_pattern, colnames(data), value = TRUE)
  
  problems <- character()
  if (length(missing_required) > 0) {
    problems <- c(problems, paste("Missing required columns:", paste(missing_required, collapse = ", ")))
  }
  if (length(isotopologues) == 0) {
    problems <- c(problems, "No isotopologue columns found (expected pattern: A+0, A+1, M+0, etc.)")
  }
  if (length(problems) > 0) {
    showNotification(paste(problems, collapse = "; "), type = "error")
    return(FALSE)
  }
  
  return(TRUE)
}

observeEvent(input$is_tracer_data, {
  if(input$is_tracer_data) {

    valid <- validate_tracer_data(rv$data[[rv$activeFile]])
    if (!valid) {  return()  }

    mfa$raw <- rv$data[[rv$activeFile]]
    metadata$metabolites <- unique(mfa$raw$Analyte)
    metadata$isotopologues <- grep("^A\\+", colnames(mfa$raw), value = TRUE) #TODO: read also M+X, other patterns?
    metadata$samples <- unique(mfa$raw$Analysis)
    


    update_choices <- function(ids, choices, type = c("select", "picker")) {
      type <- match.arg(type)
      for (id in ids) {
        tryCatch({
          if (type == "select") {
            updateSelectInput(session, id, choices = choices)
          } else {
            updatePickerInput(session, id, choices = choices)
          }
        }, error = function(e) NULL) # ignore missing UI elements
      }
    }
    update_choices(c("fc_metabolite", "ip_metabolite", "it_metabolite", "gt_metabolite", "metabolite_time_table", "mp_metabolite"), metadata$metabolites, type = "select")
    update_choices("ip_sample", metadata$samples, type = "select")
    update_choices("it_isotopologues", metadata$isotopologues, type = "picker")


    #TODO: deal with threshold since it has to be applied before normalization


    # Normalization
    numerical_cols <- mfa$raw[, metadata$isotopologues, drop = FALSE]

    normalized_sum <- numerical_cols / rowSums(numerical_cols, na.rm = TRUE)
    mfa$normalized_sum <- cbind(mfa$raw[, c("Analyte", "Analysis")], normalized_sum)

    # Normalize each row by the first isotopologue (typically A+0)
    normalized_zero <- sweep(
      numerical_cols,
      MARGIN = 1,
      STATS = mfa$raw[[metadata$isotopologues[1]]],
      FUN = "/"
    )
    print(head(normalized_zero))
    mfa$normalized_zero <- cbind(mfa$raw[, c("Analyte", "Analysis")], normalized_zero)
  }
})


observeEvent(input$inputTracerSequence, {
    valid <- validate_tracer_data(rv$data[[rv$activeFile]])
    if (!valid) {  return()  }
  
    mfa$raw <- rv$data[[rv$activeFile]]
    metadata$metabolites <- unique(mfa$raw$Analyte)
    metadata$isotopologues <- grep("^A\\+", colnames(mfa$raw), value = TRUE) #TODO: read also M+X, other patterns?
    metadata$samples <- unique(mfa$raw$Analysis)
    

    update_choices <- function(ids, choices, type = c("select", "picker")) {
      type <- match.arg(type)
      for (id in ids) {
        tryCatch({
          if (type == "select") {
            updateSelectInput(session, id, choices = choices)
          } else {
            updatePickerInput(session, id, choices = choices)
          }
        }, error = function(e) NULL) # ignore missing UI elements
      }
    }
    update_choices(c("fc_metabolite", "ip_metabolite", "it_metabolite", "gt_metabolite", "metabolite_time_table", "mp_metabolite"), metadata$metabolites, type = "select")
    update_choices("ip_sample", metadata$samples, type = "select")
    update_choices("it_isotopologues", metadata$isotopologues, type = "picker")


    #TODO: deal with threshold since it has to be applied before normalization


    # Normalization
    numerical_cols <- mfa$raw[, metadata$isotopologues, drop = FALSE]
    normalized_sum <- numerical_cols / rowSums(numerical_cols, na.rm = TRUE)
    mfa$normalized_sum <- cbind(mfa$raw[, c("Analyte", "Analysis")], normalized_sum)

    # Normalize each row by the first isotopologue (typically A+0)
    normalized_zero <- sweep(
      numerical_cols,
      MARGIN = 1,
      STATS = mfa$raw[[metadata$isotopologues[1]]],
      FUN = "/"
    )
    mfa$normalized_zero <- cbind(mfa$raw[, c("Analyte", "Analysis")], normalized_zero)

  sequence <- tryCatch(
    read.csv(input$inputTracerSequence$datapath, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(sequence)) {
    showNotification("Failed to read tracer sequence file. Please upload a valid CSV/TXT.", type = "error")
    return()
  }

  # normalize column names (case-insensitive) and require sample/group/time
  colmap <- tolower(colnames(sequence))
  required <- c("sample", "group", "time")
  missing_cols <- setdiff(required, colmap)
  if (length(missing_cols) > 0) {
    showNotification(paste0("Tracer sequence is missing required columns: ", paste(missing_cols, collapse = ", ")), type = "error")
    return()
  }

  # rename columns to a consistent casing so server code can access them reliably
  names(sequence)[which(colmap == "sample")] <- "sample"
  names(sequence)[which(colmap == "group")] <- "group"
  names(sequence)[which(colmap == "time")] <- "time"
  sequence <- sequence[, c("sample", "group", "time")]

  # check that samples in sequence are present in raw_data (also the other way around?)
  if (!all(sequence$sample %in% mfa$raw$Analysis)) {
    showNotification("Some samples in the tracer sequence are not present in the raw data.", type = "error")
    return()
  }
  if (!all(mfa$raw$Analysis %in% sequence$sample)) {
    showNotification("Some samples in the raw data are not present in the tracer sequence.", type = "error")
    return()
  }

  # update reactive values and UI
  mfa$sequence <- sequence
  metadata$groups <- unique(sequence$group)
  metadata$time_points <- sort(unique(sequence$time))
  metadata$group_time <- unique(paste(sequence$group, sequence$time, sep = "_")) #TODO should be sorted by group then time
  mfa$sequence$group_time <- paste(sequence$group, sequence$time, sep = "_")

  updateSelectInput(session, "fc_group", choices = metadata$groups)
  updateSelectInput(session, "it_group", choices = metadata$groups)
  updateSelectInput(session, "gt_group", choices = metadata$groups)
  updateSelectInput(session, "mp_group", choices = metadata$groups)
  updateSelectInput(session, "time_point", choices = metadata$time_points)
  updatePickerInput(session, "it_group_time", choices =  metadata$group_time)

  output$tracer_sequence <- renderDT({
    mfa$sequence
  })

  mfa$normalized_long_format <- format4ggplot(mfa$normalized_sum, mfa$sequence, metadata$isotopologues)

  mfa$normalized_long_format_ref <- format4ggplot(mfa$normalized_zero, mfa$sequence, metadata$isotopologues)
  
  mfa$raw_long_format <- format4ggplot(mfa$raw, mfa$sequence, metadata$isotopologues)

  #output$ggplotdata <- renderDT({
   # mfa$normalized_long_format
  #})
  output$tracer_sequence <- renderDT({ mfa$sequence })
  output$ggplotdata <- renderDT({ mfa$normalized_long_format })

})


### OVERVIEW PANEL ###

observeEvent(input$update_threshold, {
  req(input$intensity_threshold, mfa$raw)
  threshold <- as.numeric(input$intensity_threshold)
  if (!is.null(mfa$raw)) {
    raw_rows <- paste(mfa$raw$Analyte, mfa$raw$Analysis, sep = ": ")

    filtered_data <- mfa$raw[rowSums(mfa$raw[, metadata$isotopologues], na.rm = TRUE) >= threshold, ]

    output$threshold_warning <- renderUI({
      if (nrow(filtered_data) < nrow(mfa$raw)) {
        # Check which samples were removed for which metabolites
        filtered_rows <- paste(filtered_data$Analyte, filtered_data$Analysis, sep = ": ")
        removed_samples <- setdiff(raw_rows, filtered_rows)
        div(
          style = "color: orange;",
          tags$strong("Removed rows:"),
          tags$ul(
            lapply(removed_samples, function(x) tags$li(x))
          )
        )
      } else {
        div(
          style = "color: green;",
          "Intensity threshold applied successfully. No data points were removed."
        )
      }
    })
    output$tracer_table <- renderDT({
      filtered_data
    })

    threshold$data <- filtered_data

    output$tracer_table <-renderDT({
      filtered_data
    })
  }
})

#TODO Update normalized data when threshold changes


#TODO Update ggplot data when norm changes (also when sequence is uploaded)
observeEvent(c(mfa$normalized_long_format, mfa$sequence), {
  req(mfa$normalized_long_format, mfa$sequence)

  output$ggplotdata <- renderDT({
    mfa$normalized_long_format
  })
})


### Fractional Contribution ###

output$fc_plot <- renderPlotly({
  req(mfa$raw)
  req(mfa$sequence)
  req(input$fc_metabolite)
  req(input$fc_group)
  req(input$fc_plot_type)

  plot_settings$metabolite_group <- input$fc_metabolite
  plot_settings$group <- input$fc_group
  plot_settings$plot_type <- input$fc_plot_type

  # Select samples in group
  group_samples <- mfa$sequence[mfa$sequence[,'group'] %in% plot_settings$group, 'sample']
  
  plot_settings$sample <- group_samples

  plotFractionalContribution(mfa$raw, mfa$sequence, plot_settings)

})

output$fc_table <- renderDT({
  req(mfa$raw)
  req(mfa$sequence)
  req(input$fc_metabolite)
  req(input$fc_group)
  req(input$fc_plot_type)

  plot_settings$metabolite_group <- input$fc_metabolite
  plot_settings$group <- input$fc_group
  plot_settings$plot_type <- input$fc_plot_type

  # Select samples in group
  group_samples <- mfa$sequence[mfa$sequence[,'group'] %in% plot_settings$group, 'sample']
  
  plot_settings$sample <- group_samples

  selectFCtable(mfa$raw, mfa$sequence, plot_settings)

})

### Isotopologue Profiles ###

# Normalized to row sums
output$ip_plot <- renderPlotly({
  req(mfa$normalized_long_format, input$ip_metabolite, input$ip_sample)

  plot_settings$metabolite <- input$ip_metabolite
  plot_settings$sample <- input$ip_sample

  plotIsotopologueDist(mfa$normalized_long_format, plot_settings)
})

# Normalized to A.0
output$ip_plot_A0 <- renderPlotly({
  req(mfa$normalized_long_format_ref, input$ip_metabolite, input$ip_sample)

  plot_settings$metabolite <- input$ip_metabolite
  plot_settings$sample <- input$ip_sample
  
  plotIsotopologueDist2(mfa$normalized_long_format_ref, plot_settings)
})


### Isotopologue Timecourse ###
observeEvent(input$it_metabolite, {
  req(input$it_metabolite)
  req(mfa$normalized_sum)
  req(metadata$isotopologues) 

  metabolite_data <- mfa$normalized_sum[mfa$normalized_sum$Analyte == input$it_metabolite, ]
  available_iso <- intersect(metadata$isotopologues, colnames(metabolite_data)) 
  # Remove isotologues that are all NA or zero
  available_iso <- available_iso[colSums(metabolite_data[, available_iso], na.rm = TRUE) > 0]

  updatePickerInput(session, "it_isotopologues", choices = available_iso, selected = character(0))
})

observeEvent(input$update_it_plot, {
  req(input$it_metabolite, input$it_group, input$it_data_type)
  
  # Select data based on data type
  if (input$it_data_type == "raw") {
    data_source <- mfa$raw_long_format
    source_name <- "Raw"
  } else if (input$it_data_type == "normalized_sum") {
    data_source <- mfa$normalized_long_format
    source_name <- "Normalized"
  } else {
    showNotification("Invalid data type selected.", type = "error")
    return()
  }

  # defensive checks
  if (is.null(data_source) || nrow(data_source) == 0) {
    showNotification(paste(source_name, "data not available for isotopologue timecourse."), type = "error")
    return()
  }
  if (is.null(metadata$isotopologues) || length(metadata$isotopologues) == 0) { 
    showNotification("No isotopologue columns available.", type = "error")
    return()
  }

  top5_isotopologues <- NULL
  if(input$it_show_top5) {
    tryCatch(
      {
        # Get the appropriate data source for mean calculation
        if (input$it_data_type == "raw") {
          metabolite_data <- mfa$raw[mfa$raw$Analyte == input$it_metabolite, ]
        } else {
          metabolite_data <- mfa$normalized_sum[mfa$normalized_sum$Analyte == input$it_metabolite, ]
        }
        
        available_iso <- intersect(metadata$isotopologues, colnames(metabolite_data))
        if (length(available_iso) == 0) {
          showNotification("No isotopologue columns found for selected metabolite.", type = "warning")
        } else {
          mean_abundances <- colMeans(metabolite_data[, available_iso], na.rm = TRUE)
          top5_isotopologues <- names(sort(mean_abundances, decreasing = TRUE))[seq_len(min(5, length(mean_abundances)))]
          updatePickerInput(session, "it_isotopologues", selected = top5_isotopologues)
        }
      },
      error = function(e) {
        showNotification("Error determining top 5 isotopologues: please select manually.", type = "error")
        top5_isotopologues <<- NULL
      }
    )
  }
  

  # Pick isotopologues to keep (prefer user selection, then top5, then all)
  iso_choices <- NULL
  if (!is.null(input$it_isotopologues) && length(input$it_isotopologues) > 0) {
    iso_choices <- input$it_isotopologues
  } else if (!is.null(top5_isotopologues) && length(top5_isotopologues) > 0) {
    iso_choices <- top5_isotopologues
  } else {
    iso_choices <- metadata$isotopologues 
  }
  if(length(iso_choices) == 0) {
    showNotification("No isotopologues selected or available for plotting.", type = "error")
    return()
  }

   # filter data safely
  data <- tryCatch({
    df <- data_source
    df <- df[df$Analyte == input$it_metabolite, ]
    df <- df[df$Group %in% input$it_group, ]
    df <- df[df$Isotopologue %in% iso_choices, ]
    df
  }, error = function(e) {
    showNotification(paste("Error filtering isotopologue data:", conditionMessage(e)), type = "error")
    return(NULL)
  })
  if (is.null(data) || nrow(data) == 0) {
    showNotification("No data after filtering.", type = "warning")
    return()
  }

  original_samples <- unique(data$Analysis)
  print(paste("Original samples count:", length(original_samples)))

  if(!is.null(input$it_group_time) && length(input$it_group_time) > 0 ) {
    valid_group_time <- intersect(input$it_group_time, unique(data$groupTime))
    if (length(valid_group_time) == 0) {
      showNotification("Selected grouping variables (group/time) do not match available data.", type = "error")
      return()
    }
  } else {
    showNotification("Please select at least one grouping variable (group/time).", type = "error")
    return()
  }

  data <- data %>%
    dplyr::select(Analyte, Isotopologue, Abundance, Analysis, groupTime) %>%
    dplyr::filter(!is.na(Abundance)) %>%
    dplyr::filter(groupTime %in% valid_group_time)

  final_samples <- unique(data$Analysis)
  print(paste("Final samples count after filtering group/time:", length(final_samples)))
  excluded_samples <- setdiff(original_samples, final_samples)
  print("Excluded samples due to missing group/time info:")
  print(excluded_samples)
  output$exclusion_warning <- renderUI({
    if(length(excluded_samples) > 0) {
      div(
        style = "color: red;",
        paste("Warning: The following samples were excluded due to missing group/time information:", paste(excluded_samples, collapse = ", "))
      )
    } else {
      NULL
    }
  })

  # Replicates per time point - only normalize if using normalized data
  summarized <- tryCatch({
    summarized_data <- data %>%
      
      dplyr::group_by(groupTime, Isotopologue) %>%
      dplyr::summarise(
        n_replicates = n_distinct(Analysis),
        mean_abundance = mean(Abundance, na.rm = TRUE),
        sd_abundance = sd(Abundance, na.rm = TRUE),
        se_abundance = sd_abundance / sqrt(n_replicates),
        .groups = 'drop'
      ) %>%
      dplyr::select(-n_replicates) %>%
      dplyr::group_by(groupTime) %>%
      dplyr::mutate(
        cum_mean = cumsum(mean_abundance),
        segment_bottom = cum_mean - mean_abundance,
        segment_center = segment_bottom + mean_abundance / 2,
        error_lower = segment_center - se_abundance,
        error_upper = segment_center + se_abundance
      ) %>%
      dplyr::ungroup()
    
    summarized_data
    
  }, error = function(e) {
    showNotification(paste("Error summarizing data:", conditionMessage(e)), type = "error")
    return(NULL)
  })
  
  if (is.null(summarized) || nrow(summarized) == 0) {
    showNotification("No summarized data to plot.", type = "warning")
    return()
  }
  
  output$it_table  <- renderDT({ data })
  output$it_table2 <- renderDT({ summarized })

  # Pass data type to plot settings
  settings <- list(
    metabolite_iso = input$it_metabolite, 
    group_time = input$it_group_time, 
    plot_type = input$it_plot_type,
    data_type = input$it_data_type
  )
  

  output$it_plot <- renderPlotly({
    tryCatch({
      plotStackedIsotopologues(summarized, settings)
    }, error = function(e) {
      showNotification(paste("Plot error:", conditionMessage(e)), type = "error")
      NULL
    })
  })

})


### Metabolite per group ###
output$mp_plot <- renderPlotly({
  req(mfa$normalized_long_format, mfa$sequence, input$mp_metabolite, input$mp_group, input$mp_plot_type)
  
  plot_settings$metabolite <- input$mp_metabolite
  plot_settings$group <- input$mp_group
  plot_settings$plot_type <- input$mp_plot_type
  
  plotIsotopologue(mfa$normalized_long_format, mfa$sequence, plot_settings)

})

### Group x Time ###
output$group_time_plot <- renderPlotly({
  req(mfa$normalized_long_format_ref, mfa$sequence, input$metabolite_time_table, input$time_point, input$plot_type_time)
  
  plot_settings$metabolite <- input$metabolite_time_table
  plot_settings$time_points <- input$time_point
  plot_settings$plot_type <- input$plot_type_time
  
  plotGroupTime(mfa$normalized_long_format_ref, mfa$sequence, plot_settings)

})

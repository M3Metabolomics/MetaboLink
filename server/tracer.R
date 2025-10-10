mfa <- reactiveValues(
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

      #TODO: groups, time points, fractional contribution

      updateSelectInput(session, "metabolite", choices = mfa$metabolites)
      updateSelectInput(session, "metabolite_iso", choices = mfa$metabolites)
      updateSelectInput(session, "metabolite_group", choices = mfa$metabolites)
      updateSelectInput(session, "metabolite_time_table", choices = mfa$metabolites)

      updateSelectInput(session, "sample", choices = mfa$samples)


      #TODO: normalize data

      
    }
})


observeEvent(input$inputTracerSequence, {
  tracer_sequence <- read.csv(input$inputTracerSequence$datapath, header = 1, stringsAsFactors = FALSE)
  mfa$tracer_sequence <- tracer_sequence
  mfa$groups <- unique(tracer_sequence[,'group'])
  mfa$time_points <- sort(unique(tracer_sequence[,'time']))
  mfa$group_time <- unique(paste(tracer_sequence[,'group'], tracer_sequence[,'time'], sep = "_"))


  updateSelectInput(session, "group", choices = mfa$groups)
  updateSelectInput(session, "group_iso", choices = mfa$groups)
  updateSelectInput(session, "time_point", choices = mfa$time_points)
  
  updatePickerInput(session, "group_time", choices =  mfa$group_time)

  updatePickerInput(session, "isotopologues", choices = mfa$isotopologues)

  #TODO validate sequence and only show relevant columns

  output$tracer_sequence <- renderDT({
    tracer_sequence
  })

})



### Panel 2 ###

output$tracer_plot_ref <- renderPlotly({
  req(mfa$tracer_data)
  req(input$metabolite)
  req(input$sample)

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

output$tracer_table <- renderDT({
  req(mfa$tracer_data)
  req(input$metabolite)
  req(input$sample)

  plot_settings$metabolite <- input$metabolite
  plot_settings$sample <- input$sample

  metabolite_data <- mfa$tracer_data[mfa$tracer_data$Analyte == plot_settings$metabolite & mfa$tracer_data$Analysis == plot_settings$sample, ]
  metabolite_data
})


### Panel 3 ###

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


### Panel Isotopologues ###

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
    print("Selecting top 5 isotopologues")
    # Get top 5 isotopologues by mean abundance
    metabolite_data <- mfa$tracer_data[mfa$tracer_data$Analyte == isotopologue_settings$metabolite_iso, ]
    mean_abundances <- colMeans(metabolite_data[, mfa$isotopologues], na.rm = TRUE)
    top5_isotopologues <- names(sort(mean_abundances, decreasing = TRUE))[1:5]
    updatePickerInput(session, "isotopologues", selected = top5_isotopologues)
  }

  # Select relevant data: metabolite, sample
  selected_samples <- mfa$tracer_sequence[mfa$tracer_sequence[,'group'] %in% isotopologue_settings$group_iso, c('sample', 'time')]

  selected <- mfa$tracer_data[mfa$tracer_data$Analyte == isotopologue_settings$metabolite_iso & mfa$tracer_data$Analysis %in% selected_samples$sample, ]

  # Keep only selected isotopologues
  selected <- selected[, c("Analyte", "Analysis", top5_isotopologues)]

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


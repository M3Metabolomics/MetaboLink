mfa <- reactiveValues(
  tracer_data = NULL,
  tracer_sequence = NULL,
  normalized_a0 = NULL,
  normalized_sum = NULL,
  metabolites = NULL,
  isotopologues = NULL,
  samples = NULL,
  groups = NULL
)

plot_settings <- reactiveValues(
  intensity_threshold = 500,
  metabolite = NULL,
  sample = NULL,
  metabolite_group = NULL,
  group = NULL,
  plot_type = "barplot",
  metabolite_time_table = NULL,
  group_time = NULL,
  plot_type_time = "barplot"

)
# TODO user select sample and analyte columns

observeEvent(input$is_tracer_data, {
    #TODO: fill UI options
    if(input$is_tracer_data) {
      mfa$tracer_data <- rv$data[[rv$activeFile]]
      mfa$metabolites <- unique(mfa$tracer_data$Analyte)
      mfa$isotopologues <- grep("^A\\+", colnames(mfa$tracer_data), value = TRUE)
      mfa$samples <- unique(mfa$tracer_data$Analysis)

      # groups, time points, fractional contribution

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

  print(mfa$time_points)

  updateSelectInput(session, "group", choices = mfa$groups)
  updateSelectInput(session, "group_iso", choices = mfa$groups)
  updateSelectInput(session, "time_point", choices = mfa$time_points)

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
  print(group_samples)
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
  print(group_samples)
  plot_settings$sample <- group_samples

  selectFCtable(mfa$tracer_data, mfa$tracer_sequence, plot_settings)
  
})
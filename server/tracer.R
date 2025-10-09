
observeEvent(input$tracer_data, {
    #TODO: fill UI options
})


observeEvent(input$tracer_sequence, {
  tracer_sequence <- read.csv(input$tracer_sequence$datapath, header = 1, stringsAsFactors = FALSE)

  output$tracer_sequence <- renderTable({
    tracer_sequence
  })
})
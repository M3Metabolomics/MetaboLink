tracerPanel <- fluidRow(
    hidden(
        div(id = "tracer_panel",
            tabsetPanel(
                tabPanel("Sequence",
                    box(width = NULL, fluidRow(
                        column(12, fileInput("inputTracerSequence", "Upload file (.txt or .csv)",
                                         accept = c("txt/csv", "text/comma-seperated-values, text/plain", ".csv"),
                                         width = "100%")),
                        column(12, downloadButton("download_tracer_example", "Download example tracer sequence CSV")),
                        column(12, DTOutput("tracer_sequence") %>% withSpinner(color="steelblue"))
                    ))
                ),
                tabPanel("Overview",
                    box(width = NULL, 
                        fluidRow(
                            column(12,
                                selectInput("intensity_threshold",  
                                  "Select Intensity Threshold:",
                                  choices = c(100, 500, 1000, 2000, 5000, 10000),
                                  selected = 500)
                        )),
                        fluidRow(
                            column(12, box(width = NULL, DTOutput ("tracer_table") %>% withSpinner  (color="steelblue"))),
                        ) 
                    )
                ),
                tabPanel("Pl2",
                    box(width = NULL, fluidRow(
                    column(12, 
                        selectInput("metabolite",
                                    "Select Metabolite:",
                                    choices = NULL,
                                    selected = ""),
                        selectInput("sample",
                                    "Select Sample:",
                                    choices = NULL,
                                    selected = "")
                    )),
                    fluidRow(
                        column(12, box(width = NULL, plotOutput("tracer_plot_ref") %>% withSpinner(color="steelblue"))),
                        column(12, box(width = NULL, plotOutput("tracer_plot_rowsum") %>% withSpinner(color="steelblue")))
                    ),
                    fluidRow(
                        column(12, box(width = NULL, DTOutput("tracer_table") %>% withSpinner(color="steelblue")))
                    )
                )),
                tabPanel("FC",
                    box(width = NULL, fluidRow( #TODO server side fill
                    column(12, 
                        selectInput("metabolite_group",
                                    "Select Metabolite:",
                                    choices = NULL,
                                    selected = ""),
                        selectInput("group",
                                    "Select Group:",
                                    choices = NULL,
                                    selected = ""),
                        pickerInput("plot_type",
                                    "Select Plot Type:",
                                    choices = c("Error bar plot" = "errorbar", "Bar plot" = "barplot"),
                                    selected = "barplot",
                                    multiple = FALSE)
                    )),
                    fluidRow(
                        column(12, box(width = NULL, plotOutput("fc_plot") %>% withSpinner(color="steelblue")))
                    ),
                    fluidRow(
                        column(12, box(width = NULL, DTOutput("fc_table") %>% withSpinner(color="steelblue")))
                    )
                )),
                tabPanel("Isotopologues",
                    box(width = NULL, fluidRow(
                        column(12, selectInput("metabolite_iso",
                                    "Select Metabolite:",
                                    choices = NULL,
                                    selected = ""),
                        selectInput("group_iso",
                                    "Select Group:",
                                    choices = NULL,
                                    selected = ""),
                        pickerInput("isotopologues",
                                    "Select Isotopologues:",
                                    choices = NULL,
                                    selected = NULL,
                                    multiple = TRUE,
                                    options = list(`actions-box` = TRUE)),
                        checkboxInput("show_top5", "Show top 5 isotopologues only", value = FALSE),
                        pickerInput("group_time",
                                    "Select Grouping Variable (group/time):",
                                    choices = NULL,
                                    selected = NULL,
                                    multiple = TRUE,
                                    options = list(`actions-box` = TRUE)),
                        pickerInput("data_type", 
                                    "Select Data Type", 
                                    choices = c("Raw" = "raw", "Normalized" = "normalizedRowSums"),
                                    selected = "raw",
                                    multiple = FALSE),
                        pickerInput("plot_type_isotopologues", 
                                    "Select Plot Type", 
                                    choices = c("Error bar plot" = "errorbar", "Bar Plot" = "barplot"),
                                    selected = "barplot",
                                    multiple = FALSE),
                        
                        actionButton("update_iso_plot", "Generate plot")
                        #TODO add warning for excluded samples
                    )),
                    fluidRow(
                        column(12, box(width = NULL, plotOutput("iso_plot") %>% withSpinner(color="steelblue")))
                    ),
                    fluidRow(
                        column(12, box(width = NULL, DTOutput("iso_table") %>% withSpinner(color="steelblue")))
                    )
                )),
                tabPanel("Groups",
                    box(width = NULL,
                      fluidRow(column(12,
                        selectInput("metabolite_time_table",
                                    "Select Metabolite:",
                                    choices = NULL,
                                    selected = ""),
                        selectInput("time_point",
                                    "Select Time Point:",
                                    choices = NULL,
                                    selected = ""),
                        pickerInput("plot_type_time", 
                                    "Select Plot Type", 
                                    choices = c("Error bar plot" = "errorbar", "Bar plot" = "barplot"),
                                    selected = "barplot",
                                    multiple = FALSE) 
                    )),
                    fluidRow(
                        column(12, box(width = NULL, plotOutput("group_time_plot") %>% withSpinner(color="steelblue")))
                    ),
                    fluidRow(
                        column(12, box(width = NULL, DTOutput("group_time_table") %>% withSpinner(color="steelblue")))
                    )
                )),
                tabPanel("Summary",
                    box(width = NULL,
                    #    fluidRow(
                    #      column(12,htmlOutput("title")),
                    #    ),
                    #    fluidRow(
                    #      column(6, uiOutput("info_ui")),
                    #      column(6, htmlOutput("cvinfo_ui"))
                    #    )
                    #)
                ))
            )
        )   
    )
)
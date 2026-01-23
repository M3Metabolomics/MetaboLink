tracerPanel <- fluidRow(
    hidden(
        div(id = "tracer_panel",
            fluidRow(box(width = NULL,
                column(12, h2("Summary")),
                column(12,
                       p("Please ensure that you have uploaded a valid tracer sequence file before proceeding with the analysis."))
            )),
            tabsetPanel(
                tabPanel("Overview",
                    box(width = NULL,
                        fluidRow(
                            column(12, fileInput("inputTracerSequence", "Upload file (.txt or .csv)",
                                            accept = c("txt/csv", "text/comma-separated-values, text/plain", ".csv"),
                                            width = "100%")),
                            column(12,
                                selectInput("intensity_threshold",  
                                  "Select Intensity Threshold:",
                                  choices = c(100, 500, 1000, 2000, 5000, 10000),
                                  selected = 500)
                            ),
                            column(12, 
                                actionButton("update_threshold", "Update Threshold")                
                            ),
                            column(12, 
                                uiOutput("threshold_warning")
                            )
                        ),
                        fluidRow(
                            column(12, 
                                h2("Metadata"),
                                DTOutput("tracer_sequence") %>% withSpinner(color="steelblue")
                            )
                        ),
                        fluidRow(
                            column(12, 
                                h2("Tracing Data Overview")
                            ),
                            column(12, DTOutput("tracer_table") %>% withSpinner (color="steelblue")),
                            column(12, DTOutput("ggplotdata") %>% withSpinner(color="steelblue"))
                        ) 
                    )
                ),
                tabPanel("Fractional Contribution",
                    box(width = NULL, fluidRow(
                    column(12, 
                        selectInput("fc_metabolite",
                                    "Select Metabolite:",
                                    choices = NULL,
                                    selected = ""),
                        selectInput("fc_group",
                                    "Select Group:",
                                    choices = NULL,
                                    selected = ""),
                        pickerInput("fc_plot_type",
                                    "Select Plot Type:",
                                    choices = c("Error bar plot" = "errorbar", "Bar plot" = "barplot"),
                                    selected = "barplot",
                                    multiple = FALSE)
                    )),
                    fluidRow(
                        column(12, box(width = NULL, plotlyOutput("fc_plot") %>% withSpinner(color="steelblue")))
                    ),
                    fluidRow(
                        column(12, box(width = NULL, DTOutput("fc_table") %>% withSpinner(color="steelblue")))
                    )
                )),
                tabPanel("Isotopologue Profiles",
                    box(width = NULL, fluidRow(
                    column(12, 
                        selectInput("ip_metabolite",
                                    "Select Metabolite:",
                                    choices = NULL,
                                    selected = ""),
                        selectInput("ip_sample",
                                    "Select Sample:",
                                    choices = NULL,
                                    selected = "")
                    )),
                    fluidRow(
                        column(12, 
                            box(width = NULL,
                                title = "Normalized to row sums", 
                                plotlyOutput("ip_plot") %>% withSpinner(color="steelblue")))
                    ),
                    fluidRow(
                        column(12, 
                            box(width = NULL,
                                title = "Normalized to A+0", 
                                plotlyOutput("ip_plot_A0") %>% withSpinner(color="steelblue")))
                    )
                )),
                tabPanel("Isotopologue Timecourse",
                    box(width = NULL, fluidRow(
                        column(12, selectInput("it_metabolite",
                                    "Select Metabolite:",
                                    choices = NULL,
                                    selected = ""),
                        selectInput("it_group",
                                    "Select Group:",
                                    choices = NULL,
                                    selected = ""),
                        pickerInput("it_isotopologues",
                                    "Select Isotopologues:",
                                    choices = NULL,
                                    selected = NULL,
                                    multiple = TRUE,
                                    options = list(`actions-box` = TRUE)),
                        checkboxInput("it_show_top5", "Show top 5 isotopologues only", value = FALSE),
                        pickerInput("it_group_time",
                                    "Select Grouping Variable (group/time):",
                                    choices = NULL,
                                    selected = NULL,
                                    multiple = TRUE,
                                    options = list(`actions-box` = TRUE)),
                        pickerInput("it_data_type",
                                    "Select Data Type:",
                                    choices = c("Raw intensities" = "raw", "Normalized (row sums)" = "normalized_sum"),
                                    selected = "raw",
                                    multiple = FALSE),
                        pickerInput("it_plot_type", 
                                    "Select Plot Type:", 
                                    choices = c("Error bar plot" = "errorbar", "Bar Plot" = "barplot"),
                                    selected = "barplot",
                                    multiple = FALSE),
                        
                        actionButton("update_it_plot", "Generate plot"),
                        uiOutput("exclusion_warning") #TODO
                    )),
                    fluidRow(
                        column(12, box(width = NULL, plotlyOutput("it_plot") %>% withSpinner(color="steelblue")))
                    ),
                    fluidRow(
                        column(12, box(width = NULL, DTOutput("it_table") %>% withSpinner(color="steelblue")))
                    ),
                    fluidRow(
                        column(12, box(width = NULL, DTOutput("it_table2") %>% withSpinner(color="steelblue")))
                    )
                )),
                tabPanel("Metabolite Pool",
                  box(width = NULL, 
                    fluidRow(column(12, 
                        selectInput("mp_metabolite",
                                    "Select Metabolite:",
                                     choices = NULL,
                                     selected = ""),
                        selectInput("mp_group",
                                    "Select Group:",
                                     choices = NULL,
                                     selected = ""),
                        pickerInput("mp_plot_type",
                                    "Select Plot Type:",
                                     choices = c("Error bar plot"= "errorbar", "Bar plot" = "barplot"),
                                     selected = "barplot",
                                     multiple = FALSE)
                      )),
                      fluidRow(
                        column(12, box(width = NULL, plotlyOutput("mp_plot") %>% withSpinner(color="steelblue")))
                    )
                )),
                
                tabPanel("Group x Time",
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
                        column(12, box(width = NULL, plotlyOutput("group_time_plot") %>% withSpinner(color="steelblue")))
                    ),
                    fluidRow(
                        column(12, box(width = NULL, DTOutput("group_time_table") %>% withSpinner(color="steelblue")))
                    )
                ))
            )
        )   
    )
)
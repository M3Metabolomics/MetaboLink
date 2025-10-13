shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)

  # Source all server files
  serverFiles <- list.files("server", pattern = "\\.R$", full.names = TRUE)
  lapply(serverFiles, function(f) source(f, local = TRUE))  

  # Global variables
  rv <- reactiveValues(data = list(), # List of data frames
                       sequence = list(), # List of sequences
                       activeFile = NULL, # Index of active file
                       results = list(), # List of PCA results
                       tmpData = NULL, # Temporary data
                       tmpSequence = NULL, # Temporary sequence
                       choices = NULL, # List of choices
                       drift_plot_select = 1, # Drift plot selection
                       info = vector("character"), # Vector of info
                       pca_results = list(), # List of PCA results
                       outlier_df = list(), # List of outlier data frames
                       identifier_df = list(), # List of identifier data frames
                       multipleLipidNamesDf = NULL) 
  
  userConfirmation <- reactiveVal(FALSE)
  disable("upload")
  
  rankings_merge <- data.frame(
    name = c("High", "Medium", "Low"),
    priority = c(1, 2, 3)
  )
  
  ##############
  # Load files #
  ##############
  
  massCorrection <- read.csv("./csvfiles/adducts.csv") # Import mass correction data
  refmet <- read.csv("./csvfiles/refmet.csv") # Import reference metabolite data
  query <- read.csv("./csvfiles/queried_properties.csv") # Import query data
  # Hidden features
  quotes <- readLines("./csvfiles/quotes.csv")
  
  ##########################
  # Window/panel selection #
  ##########################
  
  observeEvent(list(c(input$sequence, input$example, input$upload)), {
    windowselect("sequence")
  }, ignoreInit = T)
  observeEvent(input$explore, {
    windowselect("datatable")
  })
  observeEvent(input$export, {
    windowselect("export")
  })
  observeEvent(input$statistics_button, {
    windowselect("statistics")
  })


  initializeVariables <- function() {
  rv$results[[length(rv$results) + 1]] <- list()
}

  
  renderDownloadUI <- function(idPrefix, labelSuffix) {
    renderUI({
      lapply(seq_len(length(rv$choices)), function(x) {
        fluidRow(column(12, downloadLink(paste0(idPrefix, x), paste0(rv$choices[x], labelSuffix))))
      })
    })
  }

    updateDataAndSequence <- function(notificationMessage, newFileInput, suffix, additionalInfo = NULL) {
    if (is.null(rv$tmpData)) {
      showNotification(notificationMessage, type = "error")
    } else {
      if (newFileInput) {
        newIndex <- length(rv$data) + 1
        rv$data[[newIndex]] <- rv$tmpData
        rv$sequence[[newIndex]] <- rv$tmpSequence
        newName <- paste0(names(rv$data)[rv$activeFile], suffix)
        names(rv$data)[newIndex] <- newName
        if (!is.null(additionalInfo)) {
          rv$info[newIndex] <- paste(ifelse(is.na(rv$info[rv$activeFile]), "", rv$info[rv$activeFile]), additionalInfo, "\n")
        }
        initializeVariables()
        rv$activeFile <- newIndex
      } else {
        rv$data[[rv$activeFile]] <- rv$tmpData
        rv$sequence[[rv$activeFile]] <- rv$tmpSequence
        names(rv$data)[rv$activeFile] <- paste0(names(rv$data)[rv$activeFile], suffix)
        if (!is.null(additionalInfo)) {
          rv$info[rv$activeFile] <- paste(ifelse(is.na(rv$info[rv$activeFile]), "", rv$info[rv$activeFile]), additionalInfo, "\n")
        }
      }
      rv$choices <- paste(seq_along(rv$data), ": ", names(rv$data))
      rv$tmpData <- NULL
      rv$tmpSequence <- NULL
    }
  }
  
  updateSequence <- function(seq, data, identifier_column, fill_label = NA) {
    # Get the column names from the data
    data_colnames <- colnames(data)
    
    # Identify missing rows in the sequence (new columns in data)
    missing_seq_rows <- setdiff(data_colnames, rownames(seq))
    
    # Proceed to add new identifier rows
    if (length(missing_seq_rows) > 0) {
      num_new_rows <- length(missing_seq_rows)
      
      # Create new sequence entries for the missing columns
      new_seq_entries <- data.frame(
        labels = rep(fill_label, num_new_rows),
        batch = rep(NA, num_new_rows),
        order = rep(NA, num_new_rows),
        group = rep(NA, num_new_rows),
        time = rep(NA, num_new_rows),
        paired = rep(NA, num_new_rows),
        amount = rep(NA, num_new_rows),
        stringsAsFactors = FALSE,
        row.names = missing_seq_rows
      )
      
      # Find the position of the identifier_column in seq
      identifier_position <- which(rownames(seq) == identifier_column)
      
      if (length(identifier_position) == 0) {
        # identifier_column not found in seq
        # Split seq into upper_seq and lower_seq as per your needs
        # For demonstration, we'll consider the entire seq as upper_seq and leave lower_seq empty
        upper_seq <- seq
        lower_seq <- data.frame()
      } else {
        # Split seq into upper_seq and lower_seq
        upper_seq <- seq[1:identifier_position, , drop = FALSE]
        lower_seq <- seq[(identifier_position + 1):nrow(seq), , drop = FALSE]
      }
      
      # Combine upper_seq, new_seq_entries, and lower_seq
      seq <- rbind(upper_seq, new_seq_entries, lower_seq)
    }
    
    # Reorder sequence to match data columns
    seq <- seq[data_colnames, , drop = FALSE]
    # print(seq)
    return(seq)
  }

    createDownloadHandler <- function(type, fileExtension, dataFunc) {
    function(x) {
      output[[paste0("dwn_", type, x)]] <- downloadHandler(
        filename = function() {
          paste0(names(rv$data[x]), "_", type, fileExtension)
        },
        content = function(file) {
          dataToWrite <- dataFunc(rv, x)  # Pass rv and x to the data function
          if(fileExtension == ".csv") {
            write.csv(dataToWrite, file, row.names = FALSE)
          } else if(fileExtension == ".xlsx") {
            write_xlsx(dataToWrite, file)
          } else if(fileExtension == ".txt") {
            writeLines(dataToWrite, file)
          }
        }
      )
    }
  }


  observeEvent(input$group_name_format, {
    sequence <- rv$sequence[[rv$activeFile]]
    sequence$group <- gsub("[^[:alnum:]]", "", sequence$group)
    rv$sequence[[rv$activeFile]] <- sequence
    removeModal()
  })
  
  observe({
    if (length(input$dt_drift_panel_rows_selected) == 0 && rv$drift_plot_select == 1) {
      output$driftplotoutput1 <- renderPlot({
        NULL
      })
      output$driftplotoutput21 <- renderPlot({
        NULL
      })
    } else if (rv$drift_plot_select == 1) {
      for (i in 1:(length(input$dt_drift_panel_rows_selected) + 1)) {
        output[[paste0("driftplotoutput", i)]] <- renderPlot({
          NULL
        })
        output[[paste0("driftplotoutput2", i)]] <- renderPlot({
          NULL
        })
      }
      for (i in seq_along(input$dt_drift_panel_rows_selected)) {
        local({
          my_i <- i
          output[[paste0("driftplotoutput", my_i)]] <- renderPlot({
            driftplot(
              data = rv$data[[rv$activeFile]][input$dt_drift_panel_rows_selected[my_i], ],
              seq = rv$sequence[[rv$activeFile]]
            )
          })
        })
        if (input$drift_select != "None") {
          local({
            my_i <- i
            output[[paste0("driftplotoutput2", my_i)]] <- renderPlot({
              driftplot(
                data = rv$data[[which(rv$choices %in% input$drift_select)]][input$dt_drift_panel_rows_selected[my_i], ],
                seq = rv$sequence[[rv$activeFile]]
              )
            })
          })
        }
      }
    } else if (rv$drift_plot_select == 2) {
      output$cvscatterplot <- renderPlot({
        cvscatterplot(
          data = rv$data[[rv$activeFile]],
          data2 = rv$data[[which(rv$choices %in% input$drift_select)]],
          seq = rv$sequence[[rv$activeFile]],
          name1 = names(rv$data)[rv$activeFile],
          name2 = names(rv$data)[which(rv$choices %in% input$drift_select)]
        )
      })
    }
  })

  observe({
    if (length(input$dt_boxplot_panel_rows_selected > 0)) {
      for (i in seq_along(input$dt_boxplot_panel_rows_selected)) {
        local({
          my_i <- i
          output[[paste0("boxplotoutput", my_i)]] <- renderPlot({
            myboxplot(
              data = rv$data[[rv$activeFile]][input$dt_boxplot_panel_rows_selected[my_i], ],
              seq = rv$sequence[[rv$activeFile]],
              log = input$bloxplot_log,
              ylog = input$bloxplot_ylog
            )
          })
        })
      }
    }
  })

  ###################
  # Summary of data #
  ###################
  observe({ 
    if (!is.null(rv$activeFile)) {
      # Get the data and sequence
      seq <- rv$sequence[[rv$activeFile]]
      dat <- rv$data[[rv$activeFile]]
      blank_mv <- sum(is.na(dat[seq[, 1] %in% "Blank"])) +
        sum(dat[seq[, 1] %in% "Blank"] == 0, na.rm = TRUE)
      qc_mv <- sum(is.na(dat[seq[, 1] %in% "QC"])) +
        sum(dat[seq[, 1] %in% "QC"] == 0, na.rm = TRUE)
      sample_mv <- sum(is.na(dat[seq[, 1] %in% "Sample"])) +
        sum(dat[seq[, 1] %in% "Sample"] == 0, na.rm = TRUE)
      
      sdata <- dat[seq[, 1] %in% "Sample"]
      sclass <- seq[seq[, 1] %in% "Sample", ][, 4]
      
      
      # Get the missing values
      blank_mv <- sum(is.na(dat[seq[, "labels"] %in% "Blank"])) +
        sum(dat[seq[, "labels"] %in% "Blank"] == 0, na.rm = TRUE)
      qc_mv <- sum(is.na(dat[seq[, "labels"] %in% "QC"])) +
        sum(dat[seq[, "labels"] %in% "QC"] == 0, na.rm = TRUE)
      sample_mv <- sum(is.na(dat[seq[, "labels"] %in% "Sample"])) +
        sum(dat[seq[, "labels"] %in% "Sample"] == 0, na.rm = TRUE)
      
      # Subset the data for samples
      data_subset <- dat[seq[, "labels"] %in% "Sample"]
      sclass <- seq[seq[, "labels"] %in% "Sample", ][, "group"]
      # Calculate CV
      if (sum(seq$labels %in% "QC") > 0) {
        qccv <- paste0("CV in QC samples: ",
                       round(cvmean(dat[seq[, 1] %in% "QC"]), 2), "</br>")
      } else {
        qccv <- "No QC in dataset </br>"
      }
      # Calculate CV in groups
      if (sum(!is.na(sclass)) > 0) {
        classcv <- sapply(sort(unique(sclass)), function(x) {
          round(cvmean(sdata[sclass %in% x]), 2)
        })
        classcv <- sapply(seq_along(classcv), function(x) {
          paste0("CV in group ", sort(unique(sclass))[x], ": ", classcv[x], "</br>")
        })
      } else {
        classcv <- NULL
      }
      # Combine the text
      text <- c(qccv, classcv)
      output$title <- renderText({
        HTML("<h3>", names(rv$data)[rv$activeFile], "</h3>")
      })
      # Update the information UI
      output$info_ui <- renderUI({
        HTML(nrow(dat), " features.<br>",
             ncol(dat[seq[, 1] %in% "Sample"]), " samples.<br>",
             ncol(dat[seq[, 1] %in% "QC"]), " QC samples.<br>",
             ncol(dat[seq[, 1] %in% "Blank"]), " Blank samples.<br>", "<br>",
             sample_mv, " missing values in Samples<br>",
             qc_mv, " missing values in QC samples<br>",
             blank_mv, " missing values in Blank samples<br><br>")
      })
      # Update the CV information UI
      output$cvinfo_ui <- renderUI({
        HTML(text)
      })
      
      # Update statistics select input options
      groups <- na.omit(seq[, 'group'])
      time <- na.omit(seq[, 'time'])
      
      group_inputs <- c("group1", "group2",
                        "group1_time", "group2_time",
                        "group1_polystest", "group2_polystest")
      
      # Update select input options where groups are needed
      groups <- unique(seq$group[seq$labels == "Sample" & seq$group != ""])
      time <- unique(seq$time[seq$labels == "Sample" & seq$time != ""])
      
      group_inputs <- c("group1", "group2",
                        "group1_time", "group2_time",
                        "group1_polystest", "group2_polystest",
                        "group1_enrichment","group2_enrichment",
                        "group1_cirbar", "group2_cirbar",
                        "group1_vol", "group2_vol",
                        "group1_OR", "group2_OR")
      
      for (x in group_inputs) {
        # Default selected_value to NULL
        selected_value <- NULL
        
        # Check if input corresponds to 'group2' and groups has at least 2 elements
        if (grepl("group2", x) && length(groups) > 1) {
          selected_value <- groups[2]  # Use the second group
        }
        
        # Update selectInput with appropriate choices and default selection
        updateSelectInput(session, x, label = NULL, choices = groups, selected = selected_value)
      }
      
      # Update time select input options
      time_inputs <- c("time1_time", "time2_time",
                       "time1_polystest", "time2_polystest")
      
      for (x in time_inputs) {
        updateSelectInput(session, x, label = NULL, choices = time, selected = "")
      }
      
      features <- c("FC","log2FC","p.value","p.adj","AveExpr","t","B")
      feature_inputs <- c("feature_cirbar")
      
      for (x in feature_inputs) {
        updateSelectInput(session, x, label = NULL, choices = features, selected = "log2FC")
      }
      
      columns <- colnames(dat)
      column_inputs <- c("identifier_column_refmet",
                         "annotation_column_merge",
                         "name_column_lipids",
                         "name_column_annotate", 
                         "name_column_cirbar",
                         "group_column_cirbar",
                         "group_column_heatmap")
      
      for (x in column_inputs) {
        if (grepl("identifier_column_refmet", x)) {
          default_val <- if ("Structure" %in% columns) "Structure" else ""
          updateSelectInput(session, x, label = NULL, choices = columns, selected = default_val)
        } else if (grepl("group_column_cirbar", x)) {
          default_val <- if ("super_class" %in% columns) {
            "super_class"
          } else if ("lipid_class" %in% columns) {
            "lipid_class"
          } else {
            ""
          }
          updateSelectInput(session, x, label = NULL, choices = columns, selected = default_val)
          
        } else if (grepl("name_column_lipids", x)) {
          default_val <- if ("Original annotation" %in% columns) {
            "Original annotation"
          } else if ("Original.annotation" %in% columns) {
            "Original.annotation"
          } else if ("Name" %in% columns) {
            "Name"
          } else {
            ""
          }
          updateSelectInput(session, x, label = NULL, choices = columns, selected = default_val)
          
        } else if (grepl("name_column_annotate", x)) {
          default_val <- if ("Original annotation" %in% columns) {
            "Original annotation"
          } else if ("Original.annotation" %in% columns) {
            "Original.annotation"
          } else if ("Name" %in% columns) {
            "Name"
          } else {
            ""
          }
          updateSelectInput(session, x, label = NULL, choices = columns, selected = default_val)
          
        } else if (grepl("annotation_column_merge", x)) {
          default_val <- if ("ID level" %in% columns) {
            "ID level"
          } else {
            "Name"
          }
          updateSelectInput(session, x, label = NULL, choices = columns, selected = default_val)
          
        } else {
          default_val <- if ("Original annotation" %in% columns) {
            "Original annotation"
          } else if ("Name" %in% columns) {
            "Name"
          } else {
            ""
          }
          updateSelectInput(session, x, label = NULL, choices = columns, selected = default_val)
        }
      }
    }
  })
})

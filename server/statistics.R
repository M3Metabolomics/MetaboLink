  observeEvent(input$testType, {
    sequence <- rv$sequence[[rv$activeFile]]
    enable("selectTest")
    switch(input$testType,
           GroupsUnpaired = {
             if(!any(complete.cases(sequence[, 4]))) {
               sendSweetAlert(session, "Oops!", "Invalid test. Provide information on different groups/conditions.", type = "error")
               disable("selectTest")
             }
           },
           GroupsPaired = {
             if(!any(complete.cases(sequence[, 4]))) {
               sendSweetAlert(session, "Oops!", "Invalid test. Provide information on different groups/conditions.", type = "error")
               disable("selectTest")
             }
             if(!any(complete.cases(sequence[, 6]))) {
               sendSweetAlert(session, "Oops!", "Invalid test. Indicate paired samples.", type = "error")
               disable("selectTest")
             }
           },
           GroupsTimeUnpaired = {
             if(!any(complete.cases(sequence[, 4]))) {
               sendSweetAlert(session, "Oops!", "Invalid test. Provide information on different groups/conditions.", type = "error")
               disable("selectTest")
             }
             if(!any(complete.cases(sequence[, 5]))) {
               sendSweetAlert(session, "Oops!", "Invalid test. Indicate time points.", type = "error")
               disable("selectTest")
             }
           },
           GroupsMultipleTime = {
             if(any(complete.cases(sequence[, 5])) & any(complete.cases(sequence[, 6]))) {
               sequence <- sequence[sequence[, 1] %in% "Sample" & complete.cases(sequence[, 4]), ]
               group_time <- getGroupTime(sequence)
               unique_values <- unique(group_time)
               combinations <- combn(unique_values, 2)
               valid_combinations <- combinations[, apply(combinations, 2, function(cols) is_valid_combination(cols[1], cols[2]))]
               contrasts <- generate_contrasts(as.matrix(valid_combinations)) # matrix bc if it's only 1 combination, valid_combinations is not a matrix and generate_contrasts fails
               updateCheckboxGroupInput(session, "contrasts", choices = contrasts, selected = NULL)
             } else {
               sendSweetAlert(session, "Oops!", "Invalid test. No paired samples or time points in dataset.", type = "error")
               disable("selectTest")
             }
           },
           CompareToReference = {
             if(!any(complete.cases(sequence[, 4]))) { #TODO or only one group
               sendSweetAlert(session, "Oops!", "Invalid test. Provide information on different groups/conditions.", type = "error")
               disable("selectTest")
             } else {
               updateSelectInput(session, "referenceGroup", label = NULL, choices = na.omit(sequence[, 'group']))
             }
           },
           {
             print('default')
           }
    )
  }, ignoreInit = TRUE)
  
  observeEvent(input$selectTest, {
    data <- rv$data[[rv$activeFile]]
    sequence <- rv$sequence[[rv$activeFile]]
    dataset_name <- gsub("^\\d+ :\\s+", "", rv$choices[rv$activeFile])
    
    # print(head(data))
    # print(str(sequence))
    # print(dataset_name)
    
    switch(input$testType, 
           GroupsUnpaired = {
             if(input$group1 == input$group2) {
               sendSweetAlert(session, "Oops!", "Choose different groups to compare.", type="error")
             } else {
               results <- groupComparison(data, sequence, c(input$group1, input$group2))
               rv$results[[rv$activeFile]][[length(rv$results[[rv$activeFile]])+1]] <- results
               names(rv$results[[rv$activeFile]])[length(rv$results[[rv$activeFile]])] <- paste0(dataset_name,": ",input$group1, "_vs_", input$group2)
               # print(head(rv$results[[rv$activeFile]]))
             }
           },
           GroupsPaired = {
             if(input$group1 == input$group2) {
               sendSweetAlert(session, "Oops!", "Choose different groups to compare.", type="error")
             } else {
               results <- groupComparisonPaired(data, sequence, c(input$group1, input$group2))
               rv$results[[rv$activeFile]][[length(rv$results[[rv$activeFile]])+1]] <- results
               names(rv$results[[rv$activeFile]])[length(rv$results[[rv$activeFile]])] <- paste0(dataset_name,": ",input$group1, "_vs_", input$group2)
             }
           },
           GroupsTimeUnpaired = {
             groups <- c(input$group1_time, input$group2_time)
             times <- c(input$time1_time, input$time2_time)
             groupTime <- paste(groups, times, sep = "_")
             print(groupTime)
             results <- groupComparisonTime(data, sequence, groups, times)
             rv$results[[rv$activeFile]][[length(rv$results[[rv$activeFile]])+1]] <- results
             names(rv$results[[rv$activeFile]])[length(rv$results[[rv$activeFile]])] <- paste0(dataset_name,": ",input$group1, "_vs_", input$group2)
           },
           GroupsMultipleTime = { # multi-level in limma 
             data <- data[sequence[, 1] %in% c("Name", "Sample")]
             sequence <- sequence[sequence[, 1] %in% c("Sample"), ]
             
             group_time <- getGroupTime(sequence)
             group_time <- factor(group_time, exclude = NA)
             paired <- factor(sequence[, 'paired'],  exclude = NA)
             results <- pairedAnalysis(data, group_time, input$contrasts, paired)
             
             rv$results[[rv$activeFile]] <- results
           },
           CompareToReference = {
             data <- data[sequence[, 1] %in% c("Name", "Sample")]
             groups <- sequence[complete.cases(sequence[, 4]), 4]
             results <- referenceGroupComparison(data, as.numeric(input$referenceGroup), groups)
             rv$results[[rv$activeFile]][[length(rv$results[[rv$activeFile]])+1]] <- results
           },
           {
             print('default')
           }
    )
    # Render one table for each contrast
    output$results_ui <- renderUI({
      lapply(seq_along(rv$results[[rv$activeFile]]), function(i) {
        fluidRow(
          column(12, strong(names(rv$results[[rv$activeFile]])[i])),
          column(12, box(width = NULL, DTOutput(paste0("results", i))))
        )
      })
    })
    lapply(seq_along(rv$results[[rv$activeFile]]), function(i) {
      output[[paste0("results", i)]] <- renderDT(
        rv$results[[rv$activeFile]][[i]], options = list(scrollX = TRUE)
      )
    })
  })
  observeEvent(input$seq_table_cell_edit, {
    sequence <- rv$sequence[[rv$activeFile]]
    info <- input$seq_table_cell_edit
    str(info)
    i <- info$row
    j <- info$col
    v <- info$value
    if(j == 1) {
      sendSweetAlert(session, title = "Warning", text = "Column 'labels' cannot be edited", type = "warning")
    } else {
      sequence[i, j] <- v
      rv$sequence[[rv$activeFile]] <- sequence
    }
  })
  

    observeEvent(input$updateSequence, {
    if (!is.null(rv$activeFile) && !is.null(rv$tmpSequence)) {
      rv$sequence[[rv$activeFile]] <- rv$tmpSequence
      rv$tmpSequence <- NULL
    } else {
      showNotification("No changes to update", type = "message")
    }
  })


    observeEvent(input$cleanedLipidGroup, {
    show_modal_spinner(
      spin = "atom",
      color = "#0A4F8F",
      text = "Cleaning lipid names - might take a few seconds."
    )
    
    if(is.null(rv$activeFile)) {
      showNotification("No data", type = "error")
    } else {
      sequence <- rv$sequence[[rv$activeFile]]
      data <- rv$data[[rv$activeFile]]
      identifier <- input$name_column_lipids
      
      print(identifier)
      
      original_colnames <- colnames(data)
      message("Cleaning lipid names")
      
      # Extract and clean lipid names
      lipid_names <- data[[identifier]]
      
      # remove NA or empty lipidnames
      lipid_names <- lipid_names[!is.na(lipid_names) & lipid_names != ""]
      
      # only unique names from lipid_names is allowed
      lipid_names <- unique(lipid_names)
      
      multipleLipidNamesDf <- parseLipidNames(lipid_names)
      
      # Save the cleaned lipid names for download
      rv$multipleLipidNamesDf <- multipleLipidNamesDf
      
      print(head(multipleLipidNamesDf[,1:6]))
      
      # if the multipleLipidNamesDf first column is only NA then return NULL 
      if (all(is.na(multipleLipidNamesDf[[1]]))) {
        sendSweetAlert(session, "No lipid names found", "No lipid names found in the data. Try selecting another column. ", type = "error")
        remove_modal_spinner()
        return(NULL)
      }
      
      multipleLipidNamesDf <- multipleLipidNamesDf %>%
        select(Normalized.Name,	Original.Name, Species.Name, Extended.Species.Name)
      
      # join left the data and multipleLipidNamesDf by identifier and Original.Name
      data <- left_join(
        data,
        multipleLipidNamesDf,
        by = setNames("Original.Name", identifier)) %>%
        relocate(c(Normalized.Name, Species.Name, Extended.Species.Name), .after = identifier)
      
      # rename extended species name to lipid abbreviation
      data <- data %>%
        rename(Lipid.Abbreviation = Extended.Species.Name)
      
      cat("# row data", nrow(data), "\n")
      
      if (any(is.na(data$Normalized.Name))) {
        
        message("some normalized names is NA")
        
        # Define known abbreviations to ensure complete coverage
        abbreviations_list <- c(
          "Hex3Cer", "Hex2Cer", "HexCer", "GlcCer", "GalCer", "LacCer", "C1P", "S1P", "SPH",
          "PGM", "PIP", "CDCA", "UDCA", "HDCA",
          "FA", "MG", "DG", "TG", "BMP", "CL", "PA", "PC", "PE", "PG",
          "PI", "PS", "Cer", "SM", "ST", "SE", "FC", "CE", "CA", "CAR", "DCA",
          "LCA", "GCA", "TCA", "LPE", "LNAPE"
        )
        
        # Extract all valid lipid subclasses from refmet
        unique_abbre <- sort(unique(refmet$sub_class[grepl("^[A-Z]+$", refmet$sub_class)]))
        
        # Ensure abbreviations_list is combined with known refmet subclasses
        unique_abbre <- sort(unique(c(unique_abbre, abbreviations_list)))
        
        # Sort abbreviations by length (LPC before PC)
        unique_abbre <- unique_abbre[order(nchar(unique_abbre), decreasing = TRUE)]
        
        # Create abbreviation mapping table
        abbreviations <- unique(refmet[refmet$sub_class %in% unique_abbre, c("main_class", "sub_class")])
        
        # Add missing manual mappings
        manual_mappings <- data.frame(
          main_class = c("Fatty Acid", "Ether Phospholipids","Acylcarnitines", "Ether Glycerophosphocholines", "Ether Diacylglycerols",
                         "Ether Glycerophosphoethanolamines", "Diacylglycerols", "Phosphatidylserines", "Ether Phosphatidylinositol",
                         "Triacylglycerols", "Ether Triacylglycerols", "Cholesterol Esters","Sphingoid bases",
                         "glycosphingolipids","glycosphingolipids","glycosphingolipids", "Phosphatidylcholines","Sterols", "Phosphatidylglycerol"),
          sub_class = c("FA", "O-LPE", "CAR", "O-PC", "O-DG","O-PE", "DG",
                        "O-PS", "O-PI", "TG", "O-TG", "CE", "SPB", "HexCer", "Hex2Cer", "Hex3Cer", "LPC", "ST", "PG")
        )
        
        # Combine with abbreviations
        abbreviations <- bind_rows(abbreviations, manual_mappings)
        
        # Build regex pattern, ensuring LPC is before PC
        abbreviation_pattern <- paste0("\\b(", paste(unique_abbre, collapse = "|"), ")(\\(O-)?")
        
        # print(abbreviation_pattern)
        
        # subset data to only contain those without normalized names 
        data_sub <- data[is.na(data$Normalized.Name),]
        
        data_sub <- as.data.frame(data_sub[,all_of(identifier)])

        # remove duplicated rows
        data_sub <- as.data.frame(data_sub[!duplicated(data_sub),])
        names(data_sub) <- c(identifier)
        
        data_sub <- data_sub %>%
          mutate(lipid_abbreviation = str_extract(.data[[identifier]], abbreviation_pattern))
        
        calculate_sum_name <- function(lipid_name) {
          # Extract the head group (assumed to be the first word)
          head_group <- str_extract(lipid_name, "^[A-Za-z]+")
          
          # Define a regex pattern that captures:
          #   (1) an optional "O-" prefix,
          #   (2) the carbon count (digits),
          #   (3) the double bond count (digits),
          #   (4) an optional oxygen info starting with ";" (if present),
          #   (5) an optional extra trailing info (e.g. "-SN1")
          chain_pattern <- "(O-)?(\\d+):(\\d+)(;[^\\s\\)\\/_]+)?(-[A-Za-z0-9]+)?"
          
          chains <- str_match_all(lipid_name, chain_pattern)[[1]]
          
          # If no matches are found, return NA
          if(nrow(chains) == 0) return(NA)
          
          # According to our pattern:
          #   Column 2: optional "O-" prefix
          #   Column 3: carbons
          #   Column 4: double bonds
          #   Column 5: oxygen info (e.g. ";O2")
          #   Column 6: extra trailing info (e.g. "-SN1")
          carbons <- as.numeric(chains[,3])
          dblbonds <- as.numeric(chains[,4])
          
          total_carbons <- sum(carbons, na.rm = TRUE)
          total_db <- sum(dblbonds, na.rm = TRUE)
          
          # Special case: if there's a fatty acid specification in parentheses and head group is SM,
          # add extra 2 double bonds.
          if(grepl("\\(FA", lipid_name) && head_group == "SM") {
            total_db <- total_db + 2
          }
          
          # Collect extra information from each chain match (using the first non-empty occurrence)
          extra_info <- ""
          for(i in 1:nrow(chains)) {
            prefix    <- ifelse(is.na(chains[i,2]), "", chains[i,2])  # "O-" prefix
            oxy_info  <- ifelse(is.na(chains[i,5]), "", chains[i,5])   # oxygen info (e.g. ";O2")
            trailing  <- ifelse(is.na(chains[i,6]), "", chains[i,6])   # trailing extra info (e.g. "-SN1")
            
            combined <- paste0(prefix, oxy_info, trailing)
            if(combined != "") {
              extra_info <- combined
              break
            }
          }
          
          # Adjust placement of extra info:
          # If the extra info is exactly "O-", we want it to come before the numeric chain.
          if(extra_info == "O-") {
            summed_chain <- paste0(extra_info, total_carbons, ":", total_db)
          } else {
            summed_chain <- paste0(total_carbons, ":", total_db, extra_info)
          }
          
          # For some lipids (example: when the name contains "/" and head_group is SM),
          # you may want to force the head group to lowercase.
          if(grepl("/", lipid_name) && head_group == "SM") {
            head_group <- tolower(head_group)
          }
          
          sum_name <- paste0(head_group, " ", summed_chain)
          return(sum_name)
        }
        
        # Apply the function only when lipid_abbreviation is not NA:
        data_sub <- data_sub %>%
          mutate(sum_name = mapply(function(name, abbr) {
            if (is.na(abbr)) {
              NA_character_
            } else {
              calculate_sum_name(name)
            }
          }, .data[[identifier]], lipid_abbreviation)) %>%
          relocate(sum_name, .after = .data[[identifier]])
        
        # Fix O- notation dynamically
        data_sub$lipid_abbreviation <- ifelse(
          grepl("\\(O-", data_sub$lipid_abbreviation),
          paste0("O-", gsub("\\(O-", "", data_sub$lipid_abbreviation)),
          data_sub$lipid_abbreviation
        )
        
        # Replace NA values in abbreviations with NA
        data_sub$lipid_abbreviation[is.na(data_sub$lipid_abbreviation)] <- NA
        
        # Join with abbreviations to get lipid_class
        data_sub <- data_sub %>%
          left_join(abbreviations, by = c("lipid_abbreviation" = "sub_class"),
                    relationship = "many-to-many") %>%
          { 
            if("main_class.x" %in% names(.)) {
              rename(., main_class = main_class.x) %>% select(-main_class.y)
            } else {
              .
            }
          } %>%
          rename(lipid_class = main_class) %>%
          relocate(lipid_class, .after = lipid_abbreviation)
        
        # Replace NA in lipid_class with NA
        data_sub$lipid_class[is.na(data_sub$lipid_class)] <- NA
        
        data_sub <- data_sub %>%
          distinct_at(vars(identifier), .keep_all = TRUE)
        
        message("Data with NA in sum_name or lipid_abbreviation:")
        print(data_sub[is.na(data_sub$sum_name) | is.na(data_sub$lipid_abbreviation),])
        
      }
      
      message("Head of data file: ")
      print(head(data[,1:6]))
      
      cat("# row data", nrow(data), "\n")
      
      data <- data %>%
        left_join(data_sub,
                  by = identifier,
                  relationship = "many-to-many") %>%
        mutate(
          Normalized.Name   = coalesce(Normalized.Name, sum_name),
          Species.Name      = coalesce(Species.Name, sum_name),
          Lipid.Abbreviation = coalesce(Lipid.Abbreviation, lipid_abbreviation)
        ) %>%
        rename(Lipid.Class = lipid_class) %>% 
        relocate(Lipid.Class, .after = Lipid.Abbreviation) %>%
        select(-sum_name, -lipid_abbreviation, -Lipid.Class)
      
      cat("# row data", nrow(data), "\n")

      colnames_cleaned <- setdiff(colnames(data), original_colnames)
      
      # Debug
      print("Are data column names and sequence row names identical?")
      print(identical(colnames(rv$tmpData),
                      rownames(rv$tmpSequence)))
      
      update_modal_spinner(
        session = session,
        text = "Just a second. I'm updating your sequence and data file. Please be patient."
      )
      
      # Update sequence
      print("Updating sequence...")
      updated_seq <- updateSequence(sequence, data, colnames_cleaned, "-")
      
      # Store temporarily
      rv$tmpData <- data
      rv$tmpSequence <- updated_seq
      
      # Update sequence
      print("Updating data file...")
      updateDataAndSequence(
        notificationMessage = paste("Identifiers gathered from column:", identifier),
        newFileInput = TRUE,
        suffix = "_LipClea",
        additionalInfo = NULL
      )
      
    }
    
    remove_modal_spinner()
    sendSweetAlert(session, "Success",
                   "Lipid names have been cleaned to standardization and column with lipid group has been added.",
                   type = "success")
    message(sample(quotes, 1))
  })
  
  observeEvent(input$RemoveUnannotated, {
    show_modal_spinner(
      spin = "atom",
      color = "#0A4F8F",
      text = "Remove unannotated features."
    )
    
    if(is.null(rv$activeFile)) {
      showNotification("No data", type = "error")
    } else {
      sequence <- rv$sequence[[rv$activeFile]]
      data <- rv$data[[rv$activeFile]]
      identifier <- input$name_column_annotate
      
      print(identifier)
      
      original_colnames <- colnames(data)
      
      num_feature_before <- nrow(data)
      
      
      # remove rows with NA in identifier column
      data <- data[!is.na(data[[identifier]]),]
      # remove rows that is empty string
      data <- data[data[[identifier]] != "",]
      
      num_feature_after <- nrow(data)
      
      colnames_cleaned <- setdiff(colnames(data), original_colnames)
      
      # Debug
      print("Are data column names and sequence row names identical?")
      print(identical(colnames(rv$tmpData),
                      rownames(rv$tmpSequence)))
      
      update_modal_spinner(
        session = session,
        text = "Just a second. I'm updating your sequence and data file. Please be patient."
      )
      
      # Update sequence
      print("Updating sequence...")
      updated_seq <- updateSequence(sequence, data, colnames_cleaned, "-")
      
      # Store temporarily
      rv$tmpData <- data
      rv$tmpSequence <- updated_seq
      
      # Debug
      print("Are data column names and sequence row names identical?")
      print(identical(colnames(rv$tmpData), rownames(rv$tmpSequence)))
      
      # Update sequence
      print("Updating data file...")
      updateDataAndSequence(
        notificationMessage = paste("Unnanotated features removed from column: ", identifier),
        newFileInput = TRUE,
        suffix = "_Anno",
        additionalInfo = NULL
      )
      
    }
    
    remove_modal_spinner()
    sendSweetAlert(
      session,
      title = "Success!",
      text = paste0(
        "Unannotated features have been successfully removed.\n\n",
        "Features before: ", num_feature_before, "\n",
        "Features after: ", num_feature_after
      ),
      type = "success",
      btn_labels = "OK"
    )
    message(sample(quotes, 1))
  })
  
  observeEvent(input$editColumns, {
    showModal(
      modalDialog(
        title = "Edit columns", size = "s", easyClose = TRUE,
        footer = list(actionButton("edit_cols_confirm", "Confirm"), modalButton("Dismiss")),
        fluidRow(
          column(width = 9, h4("Column name")),
          column(width = 3, style = "text-align: left;", h4("Keep"))
        ),
        lapply(seq(ncol(rv$data[[rv$activeFile]])), function(x) {
          fluidRow(
            column(
              width = 9,
              textInput(paste0("column_edit_name", x), NULL, value = colnames(rv$data[[rv$activeFile]])[x])
            ),
            column(
              width = 3, style = "text-align: center;",
              prettyCheckbox(paste0("columns_to_keep", x), NULL, status = "info", value = T)
            ),
          )
        })
      )
    )
  })
  
  observeEvent(input$edit_cols_confirm, {
    ncol <- ncol(rv$data[[rv$activeFile]])
    column_names <- character()
    column_names <- sapply(seq(ncol), function(x) {
      input[[paste0("column_edit_name", x)]]
    })
    if(!checkDuplicates(column_names)) {
      isolate(colnames(rv$data[[rv$activeFile]]) <- column_names)
      isolate(row.names(rv$sequence[[rv$activeFile]]) <- column_names)
      keep <- sapply(seq(ncol), function(x) input[[paste0("columns_to_keep", x)]])
      rv$data[[rv$activeFile]] <- rv$data[[rv$activeFile]][, keep]
      rv$sequence[[rv$activeFile]] <- rv$sequence[[rv$activeFile]][keep, ]
    }
    removeModal()
  })
  
  observeEvent(input$editGroups, {
    sequence <- rv$sequence[[rv$activeFile]]
    unique_groups <- unique(na.omit(sequence[, 4]))
    # Generate UI elements for each group
    group_ui_elements <- lapply(seq(unique_groups), function(x) {
      group <- unique_groups[x]
      fluidRow(
        column(3, h5(group)),
        column(9,
               textInput(paste0("edit_nickname", group), NULL, value = NULL)
        ),
      )
    })
    showModal(
      modalDialog(
        title = "Edit Group Nicknames", size = "s", easyClose = TRUE,
        footer = list(actionButton("group_edit_confirm", "Confirm"), modalButton("Dismiss")),
        fluidRow(
          column(3, h4("Group")),
          column(9, h4("Nickname"))
        ), 
        do.call(tagList, group_ui_elements)
      )
    )
  })
  
  observeEvent(input$group_edit_confirm, {
    sequence <- rv$sequence[[rv$activeFile]]
    groups <- sequence[, 4]
    for (x in seq_along(groups)) {
      if (!is.na(groups[x])) {
        input_x <- input[[paste0("edit_nickname", groups[x])]]
        if (nchar(input_x) != 0 & isValidName(input_x)) {
          sequence[x, 4] <- input_x
        }
      }
    }
    rv$sequence[[rv$activeFile]] <- sequence
    removeModal()
  })
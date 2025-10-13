  # Over representation analysis (ORA) ----
  observeEvent(input$select_data_for_enrichment, {
    req(input$select_data_for_enrichment)  # Ensure a dataset is selected
    
    if (!is.null(rv$activeFile)) {
      if (input$select_data_for_enrichment == "Unsaved data") {
        data <- rv$tmpData  # Use the temporary data
      } else {
        # Get the index of the selected dataset
        sd <- which(rv$choices %in% input$select_data_for_enrichment)
        data <- rv$data[[sd]]  # Retrieve the selected dataset
      }
      
      # Extract column names from the selected dataset
      data_colnames <- colnames(data)
      
      # Define the IDs for the selectInputs you wish to update
      select_ids <- c("identifier_column", "compound_column")
      
      for (col in select_ids) {
        # Set default based on which input it is
        if (col == "identifier_column") {
          default_val <- if ("Structure" %in% data_colnames) {
            "Structure"
          } else if ("InChI" %in% data_colnames) {
            "InChI"
          } else {
            ""
          }
        } else if (col == "compound_column") {
          # Change this logic as needed. For example, default to "Compound" if present.
          default_val <- if ("Original annotation" %in% data_colnames) {
            "Original annotation"
          } else {
            data_colnames[1]  # Fallback: use the first column
          }
        }
        
        updateSelectInput(session, col, choices = data_colnames, selected = default_val)
      }
    }
    
    message(sample(quotes, 1))
  })
  observeEvent(input$run_gather_identifiers, {
    req(input$select_data_for_enrichment,
        input$identifier_column,
        input$compound_column)
    
    message("Running pathway enrichment analysis...")
    
    show_modal_spinner(
      spin = "atom",
      color = "#0A4F8F",
      text = paste("Gathering Identifiers... This may take a few minutes.")
    )
    
    # Process logic
    if (!is.null(rv$activeFile)) {
      if (input$select_data_for_enrichment == "Unsaved data") {
        data <- rv$tmpData
        seq <- rv$tmpSequence
      } else {
        sd <- which(rv$choices %in% input$select_data_for_enrichment)
        data <- rv$data[[sd]]
        seq <- rv$sequence[[sd]]
      }
      
      # Get the column name for identifiers
      identifier <- input$identifier_column
      compound <- input$compound_column
      
      message(paste0("Number of features: ", nrow(data)))
      
      subset <- data
      
      message(paste0("Number of features after filtering: ", nrow(subset)))
      
      # Set this to TRUE during development and FALSE in production
      run_development_code <- TRUE
      
      if (run_development_code) {
        # Check how many rows query has initially
        
        update_modal_spinner(
          session = session,
          text = paste("Updateing query file... Please be patient. \n (1/5)")
        )
        
        query_start <- nrow(query)
        message(paste0("Number of established features: ", query_start))
        
        desired_properties <- c(
          "MolecularFormula", "MolecularWeight", "ExactMass", "MonoisotopicMass",
          "CanonicalSMILES", "IsomericSMILES", "InChI", "InChIKey", "IUPACName"
        )
        
        all_results <- data.frame()
        chunk_size <- 5
        
        # Identify new compounds (those not already in query$Identifier) using dplyr
        new_compounds <- subset %>%
          pull(!!sym(compound)) %>%  # extract the column as a vector
          { .[!. %in% query$Identifier] } %>%  # filter out entries that are in query$Identifier
          unique()
        
        clean_names <- TRUE # Set to TRUE to clean compound names else FALSE
        if (clean_names) {
          new_compounds <- gsub(";.*", "", new_compounds)
        }
        message(paste("Number of features not in query:", length(new_compounds)))
        
        if (length(new_compounds) > 0) {
          num_rows <- length(new_compounds)
          row_indices <- seq(1, num_rows, by = chunk_size)
          print("Row indices to be processed:")
          print(row_indices)
          
          if (num_rows > 0) {
            for (start_idx in row_indices) {
              end_idx <- base::min(start_idx + chunk_size - 1, num_rows)
              print(paste("Processing features from", start_idx, "to", end_idx))
              
              compound_subset <- new_compounds[start_idx:end_idx]
              print("Current features subset:")
              print(compound_subset)
              
              props_chunk <- tryCatch({
                get_properties(
                  properties = desired_properties,
                  identifier = compound_subset,
                  namespace = "name",
                  propertyMatch = list(.ignore.case = TRUE, type = "contain")
                )
              }, error = function(e) {
                warning("Failed to get properties for compounds: ", paste(compound_subset, collapse = ", "))
                return(NULL)
              })
              
              if (is.null(props_chunk) || !inherits(props_chunk, "PubChemInstanceList")) {
                print("props_chunk is NULL or not a PubChemInstanceList, skipping this batch.")
                Sys.sleep(1) # respect rate limit
                next
              }
              
              props_retrieved <- tryCatch({
                retrieve(object = props_chunk, .to.data.frame = TRUE, .combine.all = TRUE)
              }, error = function(e) {
                warning("Failed to retrieve data for compounds: ", paste(compound_subset, collapse = ", "))
                return(data.frame())
              })
              
              if (!is.null(props_retrieved) && nrow(props_retrieved) > 0) {
                all_results <- bind_rows(all_results, props_retrieved)
              } else {
                print("No valid rows retrieved for this batch.")
              }
              
              print(Sys.time())
              # Pause to respect the 5 queries/second limit
              Sys.sleep(1)
            }
          } else {
            message("No new compounds to process; skipping property queries.")
          }
          
          
          if (nrow(all_results) > 0) {
            # use dplyr to convert columns to numeric
            all_results <- all_results %>%
              mutate_at(vars(MolecularWeight,
                             ExactMass,
                             MonoisotopicMass),
                        as.numeric)
            
            dir_path <- "~/Desktop/SDU/Cand/Master thesis /GitHub/MetaboLink-Main/csvfiles" # Adjust path as needed
            file_path <- file.path(dir_path, "queried_properties.csv")
            
            # If the file already exists, read it and combine
            if (file.exists(file_path)) {
              existing_data <- read.csv(file_path, stringsAsFactors = FALSE)
              
              # Identify new entries that are not already in existing_data by 'Identifier'
              if ("Identifier" %in% colnames(existing_data) && "Identifier" %in% colnames(all_results)) {
                new_entries <- anti_join(all_results, existing_data, by = "Identifier")
              } else {
                # If 'Identifier' is not present, just append all results
                warning("No 'Identifier' column found. Appending all results without filtering.")
                new_entries <- all_results
              }
              
              # Only bind if there are new entries
              if (nrow(new_entries) > 0) {
                combined_data <- bind_rows(existing_data, new_entries)
                # Remove duplicates based on Identifier if needed
                combined_data <- combined_data[!duplicated(combined_data$Identifier), ]
                
                print(paste("Number of rows in query after gathering identifiers:", nrow(combined_data)))
                
                write.csv(combined_data, file_path, row.names = FALSE)
                message("queried_properties.csv has been updated with new entries.")
              } else {
                message("No genuinely new entries to add. The existing queried_properties.csv remains unchanged.")
              }
              
            } else {
              # No existing file, just create a new one
              write.csv(all_results, file_path, row.names = FALSE)
              message("queried_properties.csv has been created with the new entries.")
            }
            
          } else {
            message("No results retrieved for these new compounds.")
          }
        } else {
          message("No new compounds found to query.")
        }
        
      } else {
        # If run_development_code is FALSE, do nothing here
        message("Development code is disabled. Skipping property queries.")
      }
      
      update_modal_spinner(
        session = session,
        text = paste("Updating InChI... Please be patient. \n (2/5)")
      )
      
      old_columns <- colnames(subset)
      
      # Append the information from query to subset 
      subset <- subset %>%
        left_join(
          query %>% 
            select(CID, InChI, CanonicalSMILES, IsomericSMILES, InChIKey, IUPACName),
          by = setNames("InChI", identifier),
          relationship = "many-to-many"
        )
      
      if ("InChIKey.x" %in% colnames(subset)) {
        # Merge InChIKey.x and InChIKey.y into InChIKey
        subset <- subset %>%
          mutate(InChIKey = coalesce(InChIKey.x, InChIKey.y)) %>%
          select(-c(InChIKey.x, InChIKey.y))
      }
      
      subset <- subset %>%
        relocate(all_of(c("InChIKey", "CID", "CanonicalSMILES", "IsomericSMILES", "IUPACName")), .after = identifier)
      
      # Remove duplicated columns based on InChI 
      subset <- subset %>% distinct(.keep_all = TRUE)
      
      message(paste0("Number of features after joining query: ", nrow(subset)))
      # print(head(subset))
      
      # Update InChI from local cache and online queries
      message("Running update_inchi function...")
      subset_updated <- update_inchi(subset, compound, identifier, query)
      message("Update InChI function completed.")
      
      update_modal_spinner(
        session = session,
        text = paste("Looking for CID's... Please be patient. \n (3/5)")
      )
      
      # Gather CID's after InChI update
      message("Running update_cid function...")
      data_updated <- update_cid(subset_updated, compound, identifier, query)
      # message("DATA UPDATED with CID INFO: ")
      # print(head(data_updated))
      message("Update cid function completed.")
      
      update_modal_spinner(
        session = session,
        text = paste("Adding refmet information... Please be patient. \n (4/5)")
      )
      
      refmet_updated <- refmet %>%
        filter(
          ( !is.na(pubchem_cid) & pubchem_cid %in% data_updated$CID ) |
            ( !is.na(inchi_key) & inchi_key %in% data_updated$InChIKey ))
      
      filtered_refmet <- refmet_updated %>%
        filter(pubchem_cid %in% data_updated$CID |
                 inchi_key %in% data_updated$InChIKey)
      
      filtered_refmet$pubchem_cid <- as.character(filtered_refmet$pubchem_cid)
      colnames_refmet <- colnames(filtered_refmet)
      # print(head(filtered_refmet))
      
      # For data_updated, use CID if available; otherwise, use InChIKey
      data_updated <- data_updated %>%
        mutate(join_key = coalesce(CID, InChIKey))
      
      # For filtered_refmet, use pubchem_cid if available; otherwise, use InChIKey
      filtered_refmet <- filtered_refmet %>%
        mutate(join_key = coalesce(pubchem_cid, inchi_key))
      
      # Now join by join_key. Only rows where at least one key is available and matches will join.
      final_data <- data_updated %>%
        left_join(
          filtered_refmet %>% select(all_of(colnames_refmet), join_key),
          by = "join_key",
          relationship = "many-to-many"
        )
      
      final_data <- final_data %>%
        relocate(any_of(colnames_refmet), .after = IUPACName) %>%
        select(-c(pubchem_cid, inchi_key))
      
      message(paste0("Number of features after merge with refmet: ", nrow(final_data)))
      
      update_modal_spinner(
        session = session,
        text = paste("Looking for pathways... Please be patient. \n (5/5)")
      )
      
      print(head(final_data))
      
      # add pathways to final_data 
      data_updated_list <- get_kegg_pathways(final_data)
      
      pathways_long <- data_updated_list$pathways_long
      data_updated_pathways <- data_updated_list$data_joined
      
      # data_updated_pathways <- data_updated_pathways %>%
      #   select(-join_key)
      
      # add new columns to data_updated_pathways
      new_cols <- setdiff(colnames(data_updated_pathways), old_columns)
      print(new_cols)
      
      identifier_index <- which(names(data_updated_pathways) == identifier)
      
      original_cols <- names(data_updated_pathways)[!(names(data_updated_pathways) %in% new_cols)]
      
      final_col_order <- append(original_cols, new_cols, after = identifier_index)
      
      final_data <- data_updated_pathways[, final_col_order]
      
      # print(head(final_data))
      
      common_cols <- setdiff(intersect(names(final_data), names(refmet)), c(compound, "refmet_name"))
      
      # Perform the left join using final_data[[compound]] and refmet$refmet_name.
      final_data_updated <- final_data %>%
        left_join(refmet, by = setNames("refmet_name", compound), suffix = c("", ".refmet")) %>%
        # For each overlapping column (excluding the join key columns), update with the value from refmet if available.
        mutate(across(
          .cols = all_of(common_cols),
          .fns = ~ coalesce(get(paste0(cur_column(), ".refmet")), .x)
        )) %>%
        # Remove the temporary columns from refmet (those ending in ".refmet")
        select(-ends_with(".refmet"))
      
      additional_keys <- c("lipid_name", "sum_name", "Normalized.Name", "Species.Name")  # additional keys
      
      for(key in additional_keys) {
        if(key %in% names(final_data_updated)) {
          final_data_updated <- final_data_updated %>%
            left_join(refmet, by = setNames("refmet_name", key),
                      suffix = c("", paste0(".", key))) %>%
            mutate(across(
              .cols = all_of(common_cols),
              .fns = ~ coalesce(get(paste0(cur_column(), ".", key)), .x)
            )) %>%
            select(-ends_with(paste0(".", key)))
        }
      }
      
      final_data_updated <- final_data_updated %>%
        # make the pubchem_cid as charactor
        mutate(pubchem_cid = as.character(pubchem_cid))
      
      final_data_updated <- final_data_updated %>%
        mutate(
          CID = coalesce(CID, pubchem_cid),
          InChIKey = coalesce(InChIKey, inchi_key)
        ) %>%
        select(-pubchem_cid, -inchi_key)
      
      
      message("Number of features after final merge with refmet:", nrow(final_data_updated))
      
      final_data <- final_data_updated
      
      cat("\n--- Final dataset ---\n")
      print(head(final_data))
      message(paste0("Number of features after adding pathways: ", nrow(final_data)))
      
      # make a subset of the data which contains only the identifiers
      identifier_data <- final_data[, c(compound, identifier, new_cols)]
      # remove rows where NA, "", " ", "NA" are present in the compound column
      identifier_data <- identifier_data[!is.na(identifier_data[, compound]) & identifier_data[, compound] != "" & identifier_data[, compound] != " " & identifier_data[, compound] != "NA" & identifier_data[, compound] != "N/A", ]
      
      # sum column with no NA, "", " " or "NA". This will give the number of identifiers
      # print(sapply(final_data,function(x) sum(!is.na(x) & x != "" & x != " " & x != "NA" & x != "N/A")))
      
      # Store identifiers
      rv$identifier_df <- identifier_data
      
      # Debug
      print("Are data column names and sequence row names identical?")
      print(identical(colnames(rv$tmpData), rownames(rv$tmpSequence)))
      
      # Update sequence
      print("Updating sequence...")
      updated_seq <- updateSequence(seq, final_data, identifier, "-")
      
      # Store temporarily
      rv$tmpData <- final_data
      rv$tmpSequence <- updated_seq
      
      # Debug
      print("Are data column names and sequence row names identical?")
      print(identical(colnames(rv$tmpData), rownames(rv$tmpSequence)))
      
      
      # Update dataset and notify user
      updateDataAndSequence(
        notificationMessage = paste("Identifiers gathered from column:", identifier),
        newFileInput = TRUE,
        suffix = "_id",
        additionalInfo = NULL
      )
      
      sendSweetAlert(session, "Success", "Identifiers gathered and data updated.", type = "success")
    }
    
    message(sample(quotes, 1))
    
    remove_modal_spinner()
  })
  
  # Render the identifier count table
  output$identifier_count_table <- renderDT({
    req(input$select_data_for_enrichment)
    
    # Define your desired vector of column names
    desired_cols <- c("Name", "Normalized.Name", "Species.Name" , "Lipid.Abbreviation", "Original.annotation",
                      "Original annotation", "Structure", "InChI", "InChIKey", "CID", "smiles", "CanonicalSMILES", "IsomericSMILES",
                      "IUPACName", "refmet_id", "refmet_name", "super_class", "main_class", "sub_class", "chebi_id",
                      "hmdb_id", "lipidmaps_id", "kegg_id")
    
    sd <- which(rv$choices %in% input$select_data_for_enrichment)
    data <- rv$data[[sd]]
    
    # Only include columns from your desired list that exist in the selected data
    cols_to_include <- intersect(desired_cols, colnames(data))
    req(length(cols_to_include) > 0)
    
    # Compute valid counts for each desired column in data
    valid_counts <- sapply(cols_to_include, function(col) {
      sum(!is.na(data[[col]]) & data[[col]] != "")
    })
    
    # Compute NA/empty counts for each desired column in data
    na_empty_counts <- sapply(cols_to_include, function(col) {
      sum(is.na(data[[col]]) | data[[col]] == "")
    })
    
    # Create a data frame with a row for each statistic and a column for each desired column
    df <- data.frame(
      Statistic = c("Valid Count", "NA/Empty Count"),
      matrix(c(valid_counts, na_empty_counts),
             nrow = 2, byrow = TRUE, 
             dimnames = list(NULL, cols_to_include))
    )
    
    # Render the datatable with additional options for scrolling
    DT::datatable(
      df,
      options = list(
        dom = 't',
        scrollX = TRUE,
        autoWidth = TRUE
      ),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left;',
        HTML(paste0("<b>Total numbers of features: </b>", nrow(data)))
      )
    )
  })
  
  # Create a reactive expression to store the data used in dt_table
  dataForPathEnri <- reactive({
    req(input$select_data_for_enrichment)
    
    if (!is.null(rv$activeFile)) {
      if (input$select_data_for_enrichment == "Unsaved data") {
        data <- rv$tmpData
        seq <- rv$tmpSequence
      } else {
        sd <- which(rv$choices %in% input$select_data_for_enrichment)
        data <- rv$data[[sd]]
        seq <- rv$sequence[[sd]]
      }
      
      # If the specific columns are present, select only those columns.
      if (all(c("Name", "refmet_name", "main_class") %in% colnames(data))) {
        data <- data %>% select(Name, refmet_name, main_class)
      }
      data
    }
  })
  #TODO: Might be relevant to extract single or several features to quickly identify in plots. 
  # Render the datatable 
  # output$dt_table_path <- renderDT({
  #   datatable(dataForPathEnri(),
  #             options = list(autoWidth = TRUE,
  #                            scrollY = "300px",
  #                            pageLength = 50))
  # })
  # Extract the selected rows from the dt_table:
  selectedData <- reactive({
    selectedRows <- input$dt_table_path_rows_selected
    if (length(selectedRows) > 0) {
      dataForPathEnri()[selectedRows, ]
    } else {
      NULL
    }
  })
  
  # Run enrichment analysis
  observeEvent(input$run_enrichment_analysis, {
    req(input$select_data_for_enrichment,
        input$group1_enrichment,
        input$group1_enrichment,
        input$top_x_enrich)
    
    
    # print("Selected rows:")
    # print(selectedData())
    
    # Make sure group1_enrichment and group2_enrichment is not the same 
    if (input$group1_enrichment == input$group2_enrichment) {
      sendSweetAlert(session, "Error", "Select different groups for comparison.", type = "error")
      return()
    }

    # Ensure only one of gene/module is selected
    if (input$gene_selected && input$module_selected) {
      sendSweetAlert(session, "Error", "Select only one: Pathway or Module.", type = "error")
      return()
    } else if (!input$gene_selected && !input$module_selected) {
      sendSweetAlert(session, "Error", "Please select either Pathway or Module.", type = "error")
      return()
    }
    
    sd <- which(rv$choices == input$select_data_for_enrichment)
    data <- rv$data[[sd]]
    seq <- rv$sequence[[sd]]
    
    # display warning if column "kegg_id" or "kegg id" is not present
    if (!("kegg_id" %in% colnames(data) || "kegg id" %in% colnames(data))) {
      sendSweetAlert(session, "Error",
                     "No KEGG ID column found in the dataset. 
      Make sure column named 'kegg_id' is present by running 'Gather Identifiers' ", type = "error")
      return()
    }
    
    # display warning if column "refmet_name" or "refmet name" is not present
    # if (!("refmet_name" %in% colnames(data) || "refmet name" %in% colnames(data))) {
    #   sendSweetAlert(session, "Error",
    #                  "No Reference Metabolite Name column found in the dataset. 
    #   Make sure column named 'refmet_name' is present by running 'Gather Identifiers' ", type = "error")
    #   return()
    # }
    
    # display warning if column "Original annotation" is not present 
    # if (!("Original annotation" %in% colnames(data))) {
    #   sendSweetAlert(session, "Error",
    #                  "No Original annotation column found in the dataset. 
    #   Make sure column named 'Original annotation' is present before running ORA ", type = "error")
    #   return()
    # }
    
    show_modal_spinner(
      spin = "atom",
      color = "#0A4F8F",
      text = "Running analysis"
    )
    
    # Extract user-selected values
    group_of_interest <- input$group1_enrichment
    comparison_group <- input$group2_enrichment
    top_x_features <- as.numeric(input$top_x_enrich)
    dataset_name <- names(rv$data)[sd]
    
    
    p_value_thresh <- input$p_value_threshold_enrich
    pAdjustMethod <- input$pAdjustMethod_enrich
    minGSSize <- input$minGSSize_enrich
    maxGSSize <- input$maxGSSize_enrich
    qvalueCutoff <- input$qvalueCutoff_enrich
    color_con <- input$color_con_enrich
    
    message("P-value theshold: ", p_value_thresh)
    message("pAdjustMethod: ", pAdjustMethod)
    message("minGSSize: ", minGSSize)
    message("maxGSSize: ", maxGSSize)
    message("qvalueCutoff: ", qvalueCutoff)
    message("color_con: ", color_con)
    
    message("Running pathway enrichment analysis...")
    message("Selected dataset: ", dataset_name)
    message("Group of interest: ", group_of_interest)
    message("Comparison group: ", comparison_group)
    message("Top X features: ", top_x_features)
    
    # Subset sequence and data based on "Sample" labels
    filtered_sequence <- seq[seq[, "labels"] %in% c("Sample"), ]
    # Define preference order
    prefs <- c("refmet_name", "Original.annotation", "Original annotation", "Name", "name")
    
    # Find the first one that actually exists in your data.frame
    chosen <- prefs[prefs %in% names(data)][1]
    
    filtered_data <- data[, c(chosen, "kegg_id", rownames(filtered_sequence)), drop = FALSE]
    
    filtered_data <- filtered_data %>%
      filter(!is.na(kegg_id), kegg_id != "", kegg_id != " ", kegg_id != "NA", kegg_id != "N/A")
    
    message("Filtered data:")
    print(head(filtered_data))
    message("Filtered sequence:")
    print(head(filtered_sequence))
    
    selected_labels <- as.character(filtered_data[["kegg_id"]])
    fallback <- if ("Original.annotation" %in% colnames(data)) {
      as.character(data[["Original.annotation"]])
    } else if ("refmet_name" %in% colnames(data)) {
      as.character(data[["Name"]])
    } else if ("Name" %in% colnames(data)) {
      as.character(data[["Name"]])
    } else if ("name" %in% colnames(data)) {
      as.character(data[["name"]])
    } else {
      NULL
    }
    if (is.null(fallback)) {
      showNotification("No fallback column ('Name' or 'name') available.", type = "error")
      return()
    }
    missing <- is.na(selected_labels) | selected_labels == ""
    selected_labels[missing] <- fallback[missing]
    rownames(filtered_data) <- make.unique(selected_labels)

    stat_results <- calculate_stats(filtered_data, filtered_sequence, adjust.method = "BH")
    
    message("Stat results:")
    print(head(stat_results))
    
    target_contrast   <- paste0(group_of_interest, "_vs_", comparison_group)
    reversed_contrast <- paste0(comparison_group, "_vs_", group_of_interest)
    
    message("target_contrast")
    print(target_contrast)
    message("reversed_contrast")
    print(reversed_contrast)
    
    # If  contrast is present in stat_results:
    if (target_contrast %in% stat_results$Contrast) {
      # Just subset
      stats_df <- subset(stat_results, Contrast == target_contrast)
    } else if (reversed_contrast %in% stat_results$Contrast) {
      # Subset, then flip the log2FC & FC
      stats_df <- subset(stat_results, Contrast == reversed_contrast)
      stats_df$log2FC <- -stats_df$log2FC
      stats_df$FC     <- 1 / stats_df$FC
      
      # Rename the contrast column to reflect the new direction
      stats_df$Contrast <- target_contrast
    } else {
      # No matching contrast found; handle how you like (warn user, or return empty)
      warning(
        paste0("No matching contrast found for '", numerator, " vs ", denominator, "'. ",
               "Available contrasts are: ", paste(unique(stat_results$Contrast), collapse=", "))
      )
      stats_df <- data.frame()
    }
    
    stats_df <- stats_df %>%
      rownames_to_column("kegg_id") %>%
      relocate(kegg_id, .before = "Contrast")
    
    # make the rownames as the rownames of the sub_df
    rownames(stats_df) <- stats_df$kegg_id
    stats_df$kegg_id <- gsub("^[^.]+\\.", "", stats_df$kegg_id)
    
    stats_df <- stats_df %>%
      mutate(diffexpressed = case_when(
        log2FC > 0.585 & p.adj < 0.05 ~ "Upregulated", #TODO fix such that user defines
        log2FC < 0.585 & p.adj < 0.05 ~ "Downregulated", # TODO fix such that user defines
        p.adj > 0.05 ~ "Non-Significant"
      ))
    
    # remove .X from the kegg_id
    stats_df$kegg_id <- gsub("\\.\\d+$", "", stats_df$kegg_id)
    
    message("Stats df: ")
    print(head(stats_df))
    
    # If the stats_df is only contain Non-Significant then return nothing and let the user know 
    if (all(stats_df$diffexpressed == "Non-Significant", na.rm = TRUE)) {
      sendSweetAlert(session, "Error", "No significant results found. ", type = "error")
      remove_modal_spinner()
      return()
    }
    
    stat_df_signi <- stats_df[stats_df$diffexpressed != "Non-Significant",]
    
    sig_kegg_id <- sort(stat_df_signi$kegg_id)
    print("Significant KEGG IDs:")
    print(sig_kegg_id)
    
    bg_genes <- as.data.frame(unique(data$kegg_id))
    colnames(bg_genes) <- "unique_kegg_id"
    # remove rows where NA, "", " " or "NA".
    universe <- bg_genes[!is.na(bg_genes$unique_kegg_id) &
                           bg_genes$unique_kegg_id != "" &
                           bg_genes$unique_kegg_id != " " &
                           bg_genes$unique_kegg_id != "NA" &
                           bg_genes$unique_kegg_id != "N/A", ]
    
    # print(universe)
    
    # Run enrichment analysis based on selection
    if (input$gene_selected) {
      update_modal_spinner(
        session = session,
        text = "Pathway enrichment analysis in progress... Please be patient."
      )
      
      title_name <- "Pathway Enrichment Analysis"
      print(paste0("#---", title_name, " ", dataset_name, "---#"))
      
      res_geneCentric_ORA <- enrichKEGG(
        gene = sig_kegg_id,
        organism = "cpd",
        keyType = "kegg",
        pvalueCutoff = p_value_thresh,
        pAdjustMethod = pAdjustMethod,
        universe,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        qvalueCutoff = qvalueCutoff,
        use_internal_data = FALSE
      )
      
      res_df_gene_ORA <- as.data.frame(res_geneCentric_ORA)

      # Convert rownames to a column
      res_df_gene_ORA <- res_df_gene_ORA %>%
        tibble::rownames_to_column(var = "Pathway_ID")
      
      message("Enrichment results DF1:")
      print(head(res_df_gene_ORA))
      
      # if res_df_gene_ora is empty then return nothing and let the user know
      if (nrow(res_df_gene_ORA) == 0) {
        sendSweetAlert(session, "Error", "No significant results found during ORA.", type = "error")
        remove_modal_spinner()
        return()
      }
      
      # Extract the "Upregulated"/"Downregulated" status and store in a new column
      res_df_gene_ORA <- res_df_gene_ORA %>%
        mutate(Regulation = ifelse(grepl("^Upregulated", Pathway_ID), "Upregulated", "Downregulated"),
               Pathway_ID = gsub("^(Upregulated|Downregulated)\\.", "", Pathway_ID)) %>% # Remove label from Pathway_ID
        select(-Pathway_ID) %>%
        relocate(Regulation, .after = "Count")

      # Check the updated dataframe
      message("Enrichment results DF:")
      print(head(res_df_gene_ORA))
      
      # Display message if all p.adjust values are larger than the threshold
      if (all(res_df_gene_ORA$p.adjust > p_value_thresh)) {
        sendSweetAlert(session,
                       "Warning",
                       paste("All adjusted p-values values are larger than the threshold.",
                             "The lowest adjusted p-values value is", min(res_df_gene_ORA$p.adjust),
                             ". Try setting the threshold larger."), 
                       type = "Warning")
        remove_modal_spinner()
        return()
      }
      
      if (all(res_df_gene_ORA$qvalue > qvalueCutoff)) {
        sendSweetAlert(session,
                       "Warning",
                       paste("All q-values are larger than the threshold.",
                             "The lowest q-value value is", min(res_df_gene_ORA$qvalue),
                             ". Try setting the threshold larger."), 
                       type = "Warning")
        remove_modal_spinner()
        return()
      }
      
      message("EnrichRes")
      # Define the new enrichResult object
      enrichres <- new("enrichResult",
                       result = res_df_gene_ORA,  # The data frame of enrichment results
                       organism = "cpd",  # If analyzing KEGG compounds
                       keytype = "kegg",
                       ontology = "UNKNOWN",
                       gene = universe,
                       pAdjustMethod = pAdjustMethod,
                       qvalueCutoff = qvalueCutoff,
                       readable = FALSE)
      
      message("Enrichment results after new object:")
      print(class(enrichres))
      print(enrichres)
      
      
      plots <- bar_dot_plot(enrichres, title_name, top_x_features, color_con)
      
    } else {
      
      update_modal_spinner(
        session = session,
        text = "Module enrichment analysis in progress... Please be patient."
      )
      
      title_name <- "Module Enrichment Analysis"
      print(paste0("#---", title_name, " ", dataset_name, "---#"))
      
      res_moduleCentric_ORA <- enrichMKEGG(
        gene = sig_kegg_id,
        organism = "cpd",
        keyType = "kegg",
        pvalueCutoff = p_value_thresh,
        pAdjustMethod = pAdjustMethod,
        universe,
        minGSSize = minGSSize,
        qvalueCutoff = qvalueCutoff
      )
      
      res_df_module_ORA <- as.data.frame(res_moduleCentric_ORA)
      
      print("Enrich results:")
      print(head(res_df_module_ORA))
      
      # if res_df_module_ORA is empty then return nothing and let the user know
      if (nrow(res_df_module_ORA) == 0) {
        sendSweetAlert(session, "Error", "No significant results found during ORA.", type = "error")
        remove_modal_spinner()
        return()
      }
      
      # Convert rownames to a column
      res_df_module_ORA <- res_df_module_ORA %>%
        tibble::rownames_to_column(var = "Pathway_ID")
      
      # Extract the "Upregulated"/"Downregulated" status and store in a new column
      res_df_module_ORA <- res_df_module_ORA %>%
        mutate(Regulation = ifelse(grepl("^Upregulated", Pathway_ID), "Upregulated", "Downregulated"),
               Pathway_ID = gsub("^(Upregulated|Downregulated)\\.", "", Pathway_ID)) %>% # Remove label from Pathway_ID
        select(-Pathway_ID) %>%
        relocate(Regulation, .after = "ID")
      
      # Check the updated dataframe
      message("Enrichment results DF:")
      print(head(res_df_module_ORA))
      
      # Display message if all p.adjust values are larger than the threshold
      if (all(res_df_module_ORA$p.adjust > p_value_thresh)) {
        sendSweetAlert(session,
                       "Warning",
                       paste("All adjusted p-values are larger than the threshold.",
                             "The lowest adjusted p-value value is", min(res_df_module_ORA$p.adjust),
                             ". Try setting the threshold larger."), 
                       type = "Warning")
        remove_modal_spinner()
        return()
      }
      
      if (all(res_df_module_ORA$qvalue > qvalueCutoff)) {
        sendSweetAlert(session,
                       "Warning",
                       paste("All q-values are larger than the threshold.",
                             "The lowest q-value value is", min(res_df_module_ORA$qvalue),
                             ". Try setting the threshold larger."), 
                       type = "Warning")
        remove_modal_spinner()
        return()
      }
      
      # Define the new enrichResult object
      enrichres <- new("enrichResult",
                       result = res_df_module_ORA,  # The data frame of enrichment results
                       organism = "cpd",
                       keytype = "kegg",
                       ontology = "UNKNOWN",
                       gene = universe,
                       pAdjustMethod = pAdjustMethod,
                       qvalueCutoff = qvalueCutoff,
                       readable = FALSE)
      
      message("Enrichment results after new object:")
      print(class(enrichres))
      print(enrichres)
      
      plots <- bar_dot_plot(enrichres, title_name, top_x_features, color_con)
    }
    
    output$enrichment_barplot <- renderPlot({
      plots$bar
    })
    
    output$enrichment_dotplot <- renderPlot({
      plots$dot
    })
    
    # Extract sample names belonging to the selected group
    selected_group_samples <- rownames(seq[seq$group == group_of_interest & seq$labels == "Sample",])
    
    # Define preference order
    prefs <- c("refmet_name", "Original.annotation", "Original annotation", "Name", "name")
    prefs_class <- c("super_class", "main_class", "sub_class")
    
    filtered_kegg_data <- data %>%
      filter(!is.na(kegg_id), kegg_id != "", kegg_id != " ") %>%
      distinct(kegg_id, .keep_all = TRUE) %>%
      select(
        any_of(prefs), 
        kegg_id, 
        any_of(prefs_class),
        all_of(selected_group_samples)
      )
    
    enrichres_df <- as.data.frame(enrichres)
    
    title_sig <- paste0(title_name," ", dataset_name, " ")

    update_modal_spinner(
      session = session,
      text = "Plotting... Please be patient."
    )
    
    output$enrichment_cnetplot <- renderPlot({
      tryCatch({
        NetGraphPlot(enrichres_df, filtered_kegg_data, title = title_sig,
                               layout_option = input$layout_option,
                               node_size_mult = input$node_size_mult,
                               node_text_size = input$node_text_size,
                               edge_alpha = input$edge_alpha,
                               edge_width_scale = input$edge_width_scale,
                               node_color_by = input$node_color_by) 
      }, error = function(err) {
        plot.new()
        text(0.5, 0.5, "Plot failed to generate.\n
             No significant pathways detected.\n Please check your criteria or adjust the thresholds.\n.", 
             cex = 1.2, col = "red", adj = 0.5)
      })
    })
    
    output$enrichment_table <- renderDT({
      datatable(
        enrichres_df, 
        options = list(
          pageLength = 10,
          autoWidth = FALSE,
          scrollX = TRUE
        ),
        class = "display nowrap"
      )
    })
    
    message(sample(quotes, 1))
    remove_modal_spinner()
    
  })
  observeEvent(input$addRefmet, {
    show_modal_spinner(
      spin = "atom",
      color = "#0A4F8F",
      text = "Extracting RefMet information... This may take a few minutes."
    )
    
    if(is.null(rv$activeFile)) {
      showNotification("No data", type = "error")
     } else {
      sequence <- rv$sequence[[rv$activeFile]]
      data <- rv$data[[rv$activeFile]]
      identifier <- input$identifier_column_refmet
      
      print(identifier)
      
      data_query <- data %>%
        left_join(
          query %>%
            select(
              any_of(
                c(
                  "CID", "InChI","CanonicalSMILES", "IsomericSMILES",  "InChIKey", "IUPACName"))),                                 # close select(...)
          by           = setNames("InChI", identifier),
          relationship = "many-to-many"
        )
      
      if ("InChIKey.x" %in% colnames(data_query)) {
        # Merge InChIKey.x and InChIKey.y into InChIKey
        data_query <- data_query %>%
          mutate(InChIKey = coalesce(InChIKey.x, InChIKey.y)) %>%
          select(-c(InChIKey.x, InChIKey.y))
      }
      
      data_query <- data_query %>%
        relocate(any_of(c("InChIKey", "CID", "CanonicalSMILES", "IsomericSMILES", "IUPACName")), .after = identifier)
      
      
      # Remove duplicated columns based on InChI 
      data_query <- data_query %>% distinct(.keep_all = TRUE)
      
      colnames_refmet <- setdiff(names(refmet), "inchi_key")
      
      # Step 2: Merge the result with refmet using the InChIKey
      data_final <- data_query %>%
        left_join(refmet,
                  by = c("InChIKey" = "inchi_key"),
                  relationship = "many-to-many") %>%
        relocate(any_of(colnames_refmet), .after = all_of(identifier))
      
      if (input$online_refmet) {
        print("Online RefMet")
        
        update_modal_spinner(
          session = session,
          text = "Ohh so you choose to look online? This might take an extra bit of time. Please be patient."
        )
        
        sub_data <- data_final[is.na(data_final$refmet_name),]
        sub_data <- sub_data[ !grepl("^_", sub_data[[1]]), ]
        # remove duplicated rows 
        sub_data <- sub_data[!duplicated(sub_data[,identifier]),]
        
        # subset sub_data if Original annotation exist else use Name
        if("Original annotation" %in% colnames(sub_data)) {
          sub_data <- sub_data %>% 
            select("Original annotation") %>% 
            # us gsub to remove everything after ;
            mutate(`Original annotation` = gsub(";.*", "", `Original annotation`))
        } else if ("Original.annotation" %in% colnames(sub_data)) {
            sub_data <- sub_data %>% 
              select("Original.annotation") %>% 
              # us gsub to remove everything after ;
              mutate("Original.annotation" = gsub(";.*", "", `Original.annotation`))
        } else {
          sub_data <- sub_data %>% 
            select(Name)
        }
        
        # make intermediate df for online look up 
        desired_properties <- c(
          "CanonicalSMILES","IsomericSMILES","InChI","InChIKey","IUPACName")
        
        all_results <- data.frame()
        chunk_size <- 5
        
        print(head(sub_data))
        
        if (length(sub_data) > 0) {
          num_rows <- nrow(sub_data)
          print(paste0("# of entries: ", num_rows))
          
          row_indices <- seq(1, num_rows, by = chunk_size)
          print(paste0("Row indices: ", row_indices))
          
          for (start_idx in row_indices) {
            end_idx <- base::min(start_idx + chunk_size - 1, num_rows)
            print(paste0("Index: ", start_idx, "-", end_idx))
            
            sub_data_chunk <- sub_data[start_idx:end_idx,]
            print(sub_data_chunk)
            
            props_chunk <- tryCatch({
              get_properties(
                properties = desired_properties,
                identifier = sub_data_chunk,
                namespace = "name",
                propertyMatch = list(.ignore.case = TRUE, type = "contain")
              )
            }, error = function(e) {
              warning("Failed to get properties for compounds: ", paste(sub_data_chunk, collapse = ", "))
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
              warning("Failed to retrieve data for compounds: ", paste(sub_data_chunk, collapse = ", "))
              return(data.frame())
            })
            
            if (!is.null(props_retrieved) && nrow(props_retrieved) > 0) {
              all_results <- dplyr::bind_rows(all_results, props_retrieved)
            } else {
              print("No valid rows retrieved for this batch.")
            }
            
            print(Sys.time())
            # Pause to respect the 5 queries/second limit
            Sys.sleep(1)
            
          }
        } 
        
        if(!is.null(all_results)) {
          print(head(all_results))
        }

        data_query_online <- data_final %>%
          left_join(
            all_results %>% select(InChI, CanonicalSMILES, InChIKey),
            by = setNames("InChI", identifier),
            relationship = "many-to-many"
          )
        
        # merge two columns with the same name
        data_query_online <- data_query_online %>%
          mutate(
            CanonicalSMILES = coalesce(CanonicalSMILES.x, CanonicalSMILES.y),
            InChIKey = coalesce(InChIKey.x, InChIKey.y)
          ) %>%
          select(-ends_with(".x"), -ends_with(".y")) %>%
          relocate(c("CanonicalSMILES", "InChIKey"), .after = identifier)
        
        ref_cols <- c(
          "refmet_id", "refmet_name", "super_class", "main_class", 
          "sub_class", "formula", "exactmass", "pubchem_cid", 
          "chebi_id", "hmdb_id", "lipidmaps_id", "kegg_id"
        )
        
        data_final <- data_query_online %>%
          left_join(refmet,
                    by = c("InChIKey" = "inchi_key"),
                    relationship = "many-to-many") %>%
          mutate(
            across(
              .cols = any_of(ref_cols),
              .fns = ~ coalesce(
                .data[[ paste0(cur_column(), ".x") ]],
                .data[[ paste0(cur_column(), ".y") ]]
              ),
              .names = "{.col}"
            )
          ) %>%
          # then drop all the .x and .y helper columns
          select(-ends_with(".x"), -ends_with(".y")) %>%
          relocate(any_of(c("refmet_id", "refmet_name", "super_class",
                     "main_class","sub_class", "formula",
                     "exactmass", "pubchem_cid", "chebi_id",
                     "hmdb_id", "lipidmaps_id", "kegg_id")), .after = InChIKey)
        
        data_final[data_final == ""] <- NA
        data_final <- data_final[order(data_final$super_class),]
          
      } else {
        print("Offline RefMet")
        
        update_modal_spinner(
          session = session,
          text = "Very you in a hurry? This does not take as long as online lookup. Please be patient."
        )
        
        data_final <- data_final %>%
          { 
            if ("smiles" %in% colnames(.)) {
              relocate(., CanonicalSMILES, .after = smiles)
            } else {
              relocate(., CanonicalSMILES, .after = !!sym(identifier))
            }
          } %>%
          relocate(InChIKey, .after = !!sym(identifier))
        
        data_final[data_final == ""] <- NA
        data_final <- data_final[order(data_final$super_class),]
      }
      
      data_final <- data_final %>%
        relocate(InChIKey, .after = identifier) %>%
        relocate(CID, .after = InChIKey) %>%
        relocate(IsomericSMILES, .after = CanonicalSMILES) %>%
        relocate(IUPACName, .after = IsomericSMILES)
      
      data_final$CID <- as.character(data_final$CID)
      
      common_cols <- setdiff(intersect(names(data_final), names(refmet)), c("Name", "refmet_name"))
      
      # Perform the left join using data_final[[compound]] and refmet$refmet_name.
      final_data_updated <- data_final %>%
        left_join(refmet, by = setNames("refmet_name", "Name"), suffix = c("", ".refmet")) %>%
        # For each overlapping column (excluding the join key columns), update with the value from refmet if available.
        mutate(across(
          .cols = all_of(common_cols),
          .fns = ~ coalesce(get(paste0(cur_column(), ".refmet")), .x)
        )) %>%
        # Remove the temporary columns from refmet (those ending in ".refmet")
        select(-ends_with(".refmet"))
      
      additional_keys <- c("Original annotation", "Normalized.Name", "Species.Name")  # example additional keys
      
      for(key in additional_keys) {
        if(key %in% names(final_data_updated)) {
          final_data_updated <- final_data_updated %>%
            left_join(refmet, by = setNames("refmet_name", key), suffix = c("", paste0(".", key))) %>%
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
      
      data_final <- final_data_updated
      
      
      update_modal_spinner(
        session = session,
        text = "Just a second. I'm updating your sequence and data file. Please be patient."
      )
      
      # Update sequence
      print("Updating sequence...")
      updated_seq <- updateSequence(sequence, data_final, colnames_refmet, "-")
      
      # Store temporarily
      rv$tmpData <- data_final
      rv$tmpSequence <- updated_seq
      
      # Debug
      print("Are data column names and sequence row names identical?")
      print(identical(colnames(rv$tmpData), rownames(rv$tmpSequence)))
      
      # Update sequence
      print("Updating data file...")
      updateDataAndSequence(
        notificationMessage = paste("Identifiers gathered from column:", identifier),
        newFileInput = TRUE,
        suffix = "_RefMet",
        additionalInfo = NULL
      )
      
      remove_modal_spinner()
      sendSweetAlert(session, "Success", "Identifiers gathered and data updated.", type = "success")
      message(sample(quotes, 1))
    }
  })
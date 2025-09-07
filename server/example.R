  observeEvent(input$example, {
    # Load example files
    # Negative ion mode
    data <- read.csv("./example_files/Liverfetus_lipid_neg1.csv", stringsAsFactors = FALSE)
    sequence <- read.csv("./example_files/fetus seq neg.csv", stringsAsFactors = FALSE)
    # look for column "Name" or "name" and make it the first column 
    if("Name" %in% colnames(data)) {
      data <- data[, c("Name", setdiff(colnames(data), "Name"))]
    } else if("name" %in% colnames(data)) {
      data <- data[, c("name", setdiff(colnames(data), "name"))]
    }
    
    row.names(sequence) <- sequence[, 1]
    # look for row "Name" or "name" and make it the first row
    if("Name" %in% rownames(sequence)) {
      sequence <- rbind(sequence["Name", ], sequence[-which(rownames(sequence) == "Name"), ])
    } else if("name" %in% rownames(sequence)) {
      sequence <- rbind(sequence["name", ], sequence[-which(rownames(sequence) == "name"), ])
    }
    sequence <- sequence[, -1]
    rv$sequence[[length(rv$sequence) + 1]] <- sequence
    rv$data[[length(rv$data) + 1]] <- data
    names(rv$data)[length(rv$data)] <- "Liverfetus_negative"
    initializeVariables()
    
    # Positive ion mode
    data <- read.csv("./example_files/Liverfetus_lipid_pos1.csv", stringsAsFactors = FALSE)
    sequence <- read.csv("./example_files/fetus seq pos.csv", stringsAsFactors = FALSE)
    # look for column "Name" or "name" and make it the first column 
    if("Name" %in% colnames(data)) {
      data <- data[, c("Name", setdiff(colnames(data), "Name"))]
    } else if("name" %in% colnames(data)) {
      data <- data[, c("name", setdiff(colnames(data), "name"))]
    }
    row.names(sequence) <- sequence[, 1]
    # look for row "Name" or "name" and make it the first row
    if("Name" %in% rownames(sequence)) {
      sequence <- rbind(sequence["Name", ], sequence[-which(rownames(sequence) == "Name"), ])
    } else if("name" %in% rownames(sequence)) {
      sequence <- rbind(sequence["name", ], sequence[-which(rownames(sequence) == "name"), ])
    }
    sequence <- sequence[, -1]
    rv$sequence[[length(rv$sequence) + 1]] <- sequence
    rv$data[[length(rv$data) + 1]] <- data
    names(rv$data)[length(rv$data)] <- "Liverfetus_positive"
    initializeVariables()
    rv$choices <- paste(seq_along(rv$data), ": ", names(rv$data))
    
    updateTabItems(session, "tabs", selected = "Datainput")
    show("buttons")
    updateCollapse(session, "menu", close = "Data input")
    disable("example")
  })
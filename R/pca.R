pcaplot <- function(data, sequence, islog) {
  data[data == 0] <- NA # Remove zeros
  data <- data[complete.cases(data), ] # Remove rows with NA
  sequence[!is.na(sequence)] <- sequence[!is.na(sequence)] # Remove NA
  sequence[is.na(sequence)] <- "No class" # Replace NA with "No class"
  shinyCatch({ # Catch shiny errors
      ifelse(islog, data <- t(data), data <- log(t(data))) # Log transform data
    },
    blocking_level = 'warning', # Block warnings
    shiny = FALSE # Do not run in shiny
  )
  pca_results <- prcomp(data, rank. = 2, center = T, scale = F) # PCA, change rank = NULL
  pov <- summary(pca_results)[["importance"]]["Proportion of Variance", ] # Proportion of variance
  pov <- round(pov * 100, 2) # Round to two decimals
  components <- pca_results[["x"]] # PCA components
  components <- data.frame(components) # Convert to data frame
  label <- paste0(row.names(components), ": ", sequence) # Label
  col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(sequence))) # Color
  pca <- plot_ly(components, x = ~PC1, y = ~PC2, type = "scatter",
                 mode = "markers", text = label, hoverinfo = "text",
                 color = sequence, colors = col)
  pca <- pca %>% layout(
    legend = list(title = list(text = "color")), # Legend
    plot_bgcolor = "#e5ecf6", # Plot background color
    xaxis = list( # X-axis
      title = paste0("PC1 (", pov[1], "% explained var.)"), # Title
      zerolinecolor = "#ffff", # Zero line color
      zerolinewidth = 2, # Zero line width
      gridcolor = "#ffff" # Grid color
    ),
    yaxis = list(
      title = paste0("PC2 (", pov[2], "% explained var.)"), # Title
      zerolinecolor = "#ffff", # Zero line color
      zerolinewidth = 2, # Zero line width
      gridcolor = "#ffff" # Grid color
    )
  )
  return(pca)
}
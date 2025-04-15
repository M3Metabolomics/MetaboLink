pcaplot <- function(data, group, islog) {
  # remove features with missing values
  data[data == 0] <- NA
  data <- data[complete.cases(data), ]

  group[!is.na(group)] <- group[!is.na(group)] #????
  group[is.na(group)] <- "No class"

  shinyCatch({
      # log transform data if not already log transformeds
      if (islog) {
        data <- t(data)
      } else {
        data <- log(t(data))
      }
    },
    blocking_level = 'warning',
    shiny = FALSE
  )

  prin <- prcomp(data, rank. = 2, center = T, scale = F) #TODO add these options to the UI

  # extract the proportion of variance explained by each principal component.
  pov <- summary(prin)[["importance"]]["Proportion of Variance", ]
  pov <- round(pov * 100, 2)

  # extract the principal component scores (the coordinates of the data points in the new PC space).
  components <- prin[["x"]]
  components <- data.frame(components)

  # create labels for each point by combining the row names with the group information
  label <- paste0(row.names(components), ": ", group)

  # define color palette
  col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(group)))

  # create the scatter plot using plotly
  pca <- plot_ly(components, x = ~PC1, y = ~PC2, type = "scatter", mode = "markers", text = label, hoverinfo = "text", color = group, colors = col)

  x_axis_title <- paste0("PC1 (", pov[1], "% explained var.)")
  y_axis_title <- paste0("PC2 (", pov[2], "% explained var.)")
  
  pca <- pca %>% layout(
    legend = list(title = list(text = "color")),
    plot_bgcolor = "#e5ecf6",
    xaxis = list(
      title = x_axis_title,
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor = "#ffff"
    ),
    yaxis = list(
      title = y_axis_title,
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor = "#ffff"
    )
  )

  return(pca)
}
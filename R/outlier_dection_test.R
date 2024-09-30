#---------------------------------------------------------#
#                   OUTLIER DETECTION                     # 
#---------------------------------------------------------#
#### Outlier Detection ----
#### K-means clustering ----
library(factoextra)
kmeans_clustering <- function(pca_df, k, percentile_threshold, PC_df) {
  
  # Debugging to check pca_df, k, percentile_threshold
  # cat("PCA Data Frame (first few rows):\n")
  # print(str(pca_df))  # Show first few rows of pca_df to confirm it's valid
  # cat("K (Number of Clusters):", k, "\n")
  # cat("Percentile Threshold:", percentile_threshold, "\n")
  # print(str(PC_df))
  
  # Determine the optimal number of clusters
  # wss <- fviz_nbclust(pca_df[, c("PC1", "PC2")], kmeans, method = "wss", k.max = nrow(pca_df) - 1) +
  #   scale_x_discrete(
  #     breaks = seq(1, nrow(pca_df) - 1, by = 5)) +  # Reduce number of x-axis ticks
  #   labs(title = "Optimal number of clusters (Elbow Method)") +
  #   theme_bw()
  # ggplotly(wss)
  # 
  # sil <- fviz_nbclust(pca_df[, c("PC1", "PC2")], kmeans, method = "silhouette", k.max = nrow(pca_df) - 1) +
  #   scale_x_discrete(
  #     breaks = seq(1, nrow(pca_df) - 1, by = 5)) +  # Reduce number of x-axis ticks
  #   labs(title = "Optimal number of clusters (Silhouette Method)") +
  #   theme_bw()
  # ggplotly(sil)
  # gap <- fviz_nbclust(pca_df[, c("PC1", "PC2")], kmeans, method = "gap_stat", k.max = nrow(pca_df) - 1)+
  #   scale_x_discrete(
  #     breaks = seq(1, nrow(pca_df) - 1, by = 5)) +  # Reduce number of x-axis ticks
  #   labs(title = "Optimal number of clusters (Gap Statistic Method)") +
  #   theme_bw()
  # ggplotly(gap)
  
  # Extract the optimal number of clusters from silhouette method
  # layers <- sil$layers
  # for (layer in layers) {
  #   if (!is.null(layer$data) && "xintercept" %in% names(layer$data)) {
  #     x_intercepts <- layer$data$xintercept
  #     print(x_intercepts)
  #   }
  # }
  # k <- as.integer(x_intercepts[1])
  kmeans_res <- kmeans(pca_df[, c("PC1", "PC2")], centers = k)
  cluster_labels <- kmeans_res$cluster
  centroids <- kmeans_res$centers
  
  # Calculate the Euclidean distance between each point and its assigned cluster centroid
  distances <- sqrt(rowSums((pca_df[, c("PC1", "PC2")] - centroids[cluster_labels, ])^2))
  
  # Calculate the percentile-based distance threshold
  threshold_value <- quantile(distances, percentile_threshold / 100)
  
  # Create data frame with results
  kmeans_df <- data.frame(
    PC1 = pca_df[, "PC1"],
    PC2 = pca_df[, "PC2"],
    DistanceToCentroid = distances,
    Cluster = factor(cluster_labels),
    Category = ifelse(distances > threshold_value, "Outlier", "Inlier"))
  
  print(head(kmeans_df))
  
  col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(kmeans_df$Cluster)))
  
  # Plot with hover text including rownames
  # Creating the kmeans_plot with hover text
  kmeans_plot <- ggplot(kmeans_df, aes(x = PC1, y = PC2, color = Cluster, shape = Category)) +
    geom_point(aes(text = paste("Sample: ", rownames(kmeans_df),
                                "<br>Cluster: ", Cluster,
                                "<br>Category: ", Category,
                                "<br>PC1: ", round(PC1, 2),
                                "<br>PC2: ", round(PC2, 2))), size = 2) +
    geom_point(data = as.data.frame(centroids),
               aes(text = paste("Cluster centroid: ", 1:k,
                                "<br>Centroid PC1: ", round(PC1, 2),
                                "<br>Centroid PC2: ", round(PC2, 2)),
                   x = PC1, y = PC2), color = col, size = 2, shape = 4) +
    labs(title = "K-means Clustering with Outlier Detection",
         x = paste0("PC1 (", round(PC_df[1,2], 2), "% explained var.)"),
         y = paste0("PC2 (", round(PC_df[1,2], 2), "% explained var.)")) +
    theme_bw() +
    scale_color_manual(values = col) +
    scale_shape_manual(values = c("Inlier" = 16, "Outlier" = 17))
  # Convert ggplot to plotly for interactive hover tooltips
  kmeans_plot_plotly <- ggplotly(kmeans_plot, tooltip = "text")
  
  return(list(kmeans_df = kmeans_df,
              kmeans_plotly = kmeans_plot_plotly))
}

#### 
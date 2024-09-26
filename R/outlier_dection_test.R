#---------------------------------------------------------#
#                   OUTLIER DETECTION                     # 
#---------------------------------------------------------#
#### Outlier Detection ----
#### K-means clustering ----
library(factoextra)

kmeans_clustering_with_distances <- function(pca_df, percentile_threshold) {
  # Determine the optimal number of clusters
  wss <- fviz_nbclust(pca_df[, c("PC1", "PC2")], kmeans, method = "wss", k.max = nrow(pca_df) - 1) +
    scale_x_discrete(
      breaks = seq(1, nrow(pca_df) - 1, by = 5)) +  # Reduce number of x-axis ticks
    labs(title = "Optimal number of clusters (Elbow Method)") +
    theme_bw()
  ggplotly(wss)
  
  sil <- fviz_nbclust(pca_df[, c("PC1", "PC2")], kmeans, method = "silhouette", k.max = nrow(pca_df) - 1) + 
    scale_x_discrete(
      breaks = seq(1, nrow(pca_df) - 1, by = 5)) +  # Reduce number of x-axis ticks
    labs(title = "Optimal number of clusters (Silhouette Method)") +
    theme_bw()
  ggplotly(sil)
  gap <- fviz_nbclust(pca_df[, c("PC1", "PC2")], kmeans, method = "gap_stat", k.max = nrow(pca_df) - 1)+
    scale_x_discrete(
      breaks = seq(1, nrow(pca_df) - 1, by = 5)) +  # Reduce number of x-axis ticks
    labs(title = "Optimal number of clusters (Gap Statistic Method)") +
    theme_bw()
  ggplotly(gap)
  
  # Extract the optimal number of clusters from silhouette method
  layers <- sil$layers
  for (layer in layers) {
    if (!is.null(layer$data) && "xintercept" %in% names(layer$data)) {
      x_intercepts <- layer$data$xintercept
      print(x_intercepts)
    }
  }
  
  k <- as.integer(x_intercepts[1])
  kmeans_res <- kmeans(pca_df[, 2:3], centers = k)
  cluster_labels <- kmeans_res$cluster
  centroids <- kmeans_res$centers
  
  # Calculate distances to centroids
  distances <- sqrt(rowSums((pca_df[, 2:3] - centroids[cluster_labels, ])^2))
  
  # Calculate the percentile-based distance threshold
  threshold_value <- quantile(distances, percentile_threshold / 100)
  
  # Create data frame with results
  kmeans_df <- data.frame(
    Sample = pca_df$Sample,
    PC1 = pca_df[, 2],
    PC2 = pca_df[, 3],
    DistanceToCentroid = distances,
    cluster = factor(cluster_labels),
    Category = ifelse(distances > threshold_value, "Outlier", "Inlier")
  )
  
  # Return the data frame
  return(kmeans_df)
}
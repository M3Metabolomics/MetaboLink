#---------------------------------------------------------#
#                   OUTLIER DETECTION                     # 
#---------------------------------------------------------#
# Load required libraries
library(ggplot2)
library(dbscan)
library(fpc)
library(plotly)
library(factoextra)
library(caret)
library(ggdendro)
library(reshape2)
library(car)
library(heatmaply)
library(party)
library(rpart)
library(rpart.plot)
library(ROCR)
library(randomForest)
library(Rlof)

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

#### Hierarchical clustering ----
library(ggdendro)
perform_hierarchical_clustering <- function(pca_df, sequence, method = "complete", k = 3, threshold = 3) {
  # Debugging to check and sequence
  # print(head(sequence))
  # Ensure rownames are unique
  rownames(sequence) <- ifelse(grepl("^[0-9]+$", rownames(sequence)), 
                               paste0("X", rownames(sequence)), 
                               rownames(sequence))
  
  # Debugging to check pca_df and sequence
  # print(head(pca_df))
  # print(head(sequence))
  
  # Compute distance matrix and perform hierarchical clustering
  data_dist <- dist(pca_df[, c("PC1", "PC2")])
  hc <- hclust(data_dist, method = method)
  
  # Perform the clustering with k clusters
  clusters <- cutree(hc, k = k)
  
  # Ensure clusters are in a factor that matches the 'group' levels in sequence
  cluster_factors <- factor(clusters) 
  
  # Debugging to check clusters and cluster_factors
  # cat("Clusters:\n", clusters, "\n")
  # cat("Clusters Factor:\n", cluster_factors, "\n")
  
  # Custom color palette for clusters
  col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(cluster_factors)))
  # cat("Color Palette: ", col, "\n")
  
  # Use scale_color_manual to apply the custom color palette
  dimr_hier_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster_factors)) +
    geom_point(shape = 16) +
    labs(title = paste("Hierarchical Clustering -", method, "method"), x = "PC1", y = "PC2") +
    scale_color_manual(values = col) +
    theme_bw()
  
  # Convert hclust object to dendrogram
  dend_data <- as.dendrogram(hc)
  dendro_data <- dendro_data(dend_data)
  
  # Debugging to check dendro_data
  # cat("Dendrogram Data: ", str(dendro_data), "\n")
  
  # Create hierarchical clustering dendrogram plot
  hc_plot <- ggplot(dendro_data$segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_bw() +
    labs(title = paste("Hierarchical Clustering -", method, "method"),
         x = "Data Points", 
         y = "Height (Distance)")
  
  # Match dendrogram labels to PCA sample data
  heights <- data.frame(label = dendro_data$labels$label, height = 0)
  
  # Ensure label heights are matched correctly
  for (i in seq_along(dendro_data$segments$yend)) {
    if (!is.na(dendro_data$segments$y[i])) {
      label <- dendro_data$segments$x[i]
      height <- dendro_data$segments$y[i]
      if (label <= nrow(heights)) {
        heights$height[label] <- height
      }
    }
  }
  
  # Match heights with PCA data based on sample labels
  heights_matched <- heights[match(pca_df$sample, heights$label), "height"]
  
  # Categorize samples based on threshold
  hc_outliers <- data.frame(
    Sample = pca_df$sample,
    Cluster = cluster_factors,
    Height = heights_matched,
    Category = ifelse(heights_matched > threshold, "Outlier", "Inlier")  # Use the specified threshold
  )
  
  # Return all plots and data as a list
  return(list(
    hclust_plot = ggplotly(dimr_hier_plot),
    dendrogram_plot = ggplotly(hc_plot),
    hierarchical_outliers = hc_outliers))
}


# Function for confusion matrix plot
perform_confusion_matrix <- function(pca_df, sequence, method = "complete", k = 3) {
  data_dist <- dist(pca_df[, c("PC1", "PC2")])
  hc <- hclust(data_dist, method = method)
  clusters <- cutree(hc, k = k)
  
  # Debugging to check clusters
  cat("Clusters:\n", clusters, "\n")
  
  reference <- factor(sequence$group)
  prediction <- factor(clusters, levels = levels(reference))
  # Debug to check reference and prediction
  cat("Reference:\n", reference, "\n")
  cat("Prediction:\n", prediction, "\n")
  
  cm <- confusionMatrix(prediction, reference)
  # debug to check confusion matrix
  cat("Confusion Matrix:\n", cm$table, "\n")
  
  cm_df <- as.data.frame(cm$table)
  cm_df$Reference <- factor(cm_df$Reference, levels = rev(levels(cm_df$Reference)))
  
  # Debugging to check cm_df
  cat("Confusion Matrix Data Frame:\n", cm_df, "\n")
  
  cm_plot <- ggplot(cm_df, aes(x = Prediction, y = Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "firebrick") +
    labs(x = "Prediction", y = "Reference") +
    theme_bw()
  
  ggplotly(cm_plot)
}


# Dendrogram Function
perform_dendrogram <- function(pca_df, sequence, method = "complete", threshold) {
  data_dist <- dist(pca_df[, c("PC1", "PC2")])
  hc <- hclust(data_dist, method = method)
  
  dend_data <- as.dendrogram(hc)
  dendro_data <- dendro_data(dend_data)
  
  # Highlight segments above the threshold
  dendro_data$segments$highlight <- dendro_data$segments$y > threshold | dendro_data$segments$yend > threshold
  
  dend_plot <- ggplot() + 
    geom_segment(data = dendro_data$segments, aes(x = x, y = y, xend = xend, yend = yend, color = highlight)) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +  # Red color for outliers
    geom_text(data = dendro_data$labels, aes(x = x, y = y, label = label), hjust = 1, angle = 45, size = 3) +
    scale_x_continuous(breaks = NULL) + # Remove numeric x-axis values
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_blank()) +
    labs(title = paste("Dendrogram -", method, "method"), 
         x = "Data Points", 
         y = "Height (Distance)")
  
  ggplotly(dend_plot)
}

#### DBSCAN clustering ----

perform_kNN_dist_plot <- function(PCA.df, k) {
  # Ensure k is numeric and within valid range
  if (!is.numeric(k) || k <= 0 || k > nrow(PCA.df)) {
    stop("Invalid k value for kNN distance plot")
  }
  
  # Debugging information
  print(paste("k value:", k))
  print(paste("Dimensions of PCA.df:", dim(PCA.df)))
  
  # Compute kNN distances
  kNN_dist <- dbscan::kNNdist(PCA.df[, c("PC1","PC2")], k = k)
  
  # Debugging information
  print(head(kNN_dist))
  
  # Create a data frame for plotting
  kNN_dist_df <- data.frame(Distance = sort(kNN_dist))
  kNN_dist_df$Index <- seq_len(nrow(kNN_dist_df))
  
  # Manually create the kNN distance plot using ggplot2
  kNN_plot <- ggplot(kNN_dist_df, aes(x = Index, y = Distance)) +
    geom_line() +
    geom_hline(yintercept = mean(kNN_dist), color = "red", linetype = "dashed") +
    labs(title = paste("kNN Distance Plot for k =", k), x = "Points sorted by distance", y = "kNN Distance") +
    theme_bw()
  
  return(ggplotly(kNN_plot))
}

# DBSCAN Function
perform_dbscan_clustering <- function(PCA.df, eps, minPts) {
  # Ensure eps and minPts are numeric
  if (!is.numeric(eps) || eps <= 0) {
    stop("Invalid eps value for DBSCAN")
  }
  if (!is.numeric(minPts) || minPts <= 0) {
    stop("Invalid minPts value for DBSCAN")
  }
  
  dbscan_res <- dbscan::dbscan(PCA.df[, c("PC1","PC2")], eps = eps, minPts = minPts)
  
  plot_data <- data.frame(Sample = rownames(PCA.df), 
                          PC1 = PCA.df$PC1,
                          PC2 = PCA.df$PC2,
                          Cluster = factor(dbscan_res$cluster),
                          Outlier = ifelse(dbscan_res$cluster == 0, "Outlier", "Inlier"))
  
  DBSCAN_plot <- ggplot(plot_data, aes(text = Sample, x = PC1, y = PC2, color = Cluster)) +
    geom_point() +
    labs(title = paste("DBSCAN Clustering - eps:", eps, "- minPts:", minPts), x = "PC1", y = "PC2", color = "Cluster") +
    theme_bw()
  
  list(
    dbscan_plot = ggplotly(DBSCAN_plot),
    dbscan_outliers = plot_data
  )
}
#### HDBSCAN clustering ----
library(dbscan)
perform_hdbscan_clustering <- function(PCA.df, minPts) {
  # Ensure minPts is numeric
  if (!is.numeric(minPts) || minPts <= 0) {
    stop("Invalid minPts value for HDBSCAN")
  }
  
  res_hdbscan <- hdbscan(dist(PCA.df[, c("PC1","PC2")], method = "euclidean"), minPts = minPts)
  
  plot_data <- data.frame(Sample = rownames(PCA.df), 
                          PC1 = PCA.df$PC1,
                          PC2 = PCA.df$PC2,
                          Cluster = factor(res_hdbscan$cluster),
                          OutlierScore = res_hdbscan$outlier_scores)
  
  HDBSCAN_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Cluster, text = Sample, outlierscore = OutlierScore)) +
    geom_point(pch = 19) +
    labs(title = paste("HDBSCAN Clustering - MinPts:", minPts), x = "PC1", y = "PC2", color = "Cluster", size = "Outlier Score") +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 4)), size = guide_legend(override.aes = list(alpha = 0.5)))
  
  list(
    hdbscan_plot = ggplotly(HDBSCAN_plot),
    hdbscan_outliers = plot_data
  )
}

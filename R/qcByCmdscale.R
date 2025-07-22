qcByCmdscale <- function(CrcBiomeScreenObject,
                         TaskName = NULL,
                         normalize_method = NULL,
                         threshold_sd = 1) {
  # Extract normalized data and sample IDs
  study_data <- CrcBiomeScreenObject$NormalizedData
  sampleID <- rownames(CrcBiomeScreenObject$SampleData)

  # Step 1: Compute Bray-Curtis distance matrix and apply classical MDS
  dist_data <- dist(study_data)
  mds_coords <- cmdscale(dist_data, k = 2)
  rownames(mds_coords) <- sampleID

  x <- mds_coords[, 1]
  y <- -mds_coords[, 2] # Flip axis 2 for better visualization

  # Step 2: Compute Euclidean distance of each sample to the centroid
  center <- c(mean(x), mean(y))
  distances <- sqrt((x - center[1])^2 + (y - center[2])^2)

  # Step 3: Identify outliers based on mean + n*SD rule
  threshold <- mean(distances) + threshold_sd * sd(distances)
  is_outlier <- distances > threshold
  outliers <- sampleID[is_outlier]

  # Print outlier sample IDs
  message("Outlier sample IDs detected:")
  print(outliers)

  # Ensure IsOutlier is a factor for color mapping
  plot_df <- data.frame(
    SampleID = sampleID,
    Dim1 = x,
    Dim2 = y,
    IsOutlier = factor(is_outlier, levels = c(FALSE, TRUE))
  )
  plot_df$IsOutlier <- factor(plot_df$IsOutlier, levels = c(FALSE, TRUE))
  # Output file
  pdf_name <- paste0("cmdscale_", TaskName, "_", normalize_method, ".pdf")

  # Create plot
  p <- ggplot(plot_df, aes(x = Dim1, y = Dim2)) +
    geom_point(color = "grey50", size = 2) + # 背景点为灰色
    geom_point(
      data = subset(plot_df, IsOutlier == TRUE),
      aes(x = Dim1, y = Dim2),
      color = "red", size = 2, shape = 1, stroke = 1.2
    ) + # 离群点为红色
    geom_text(aes(label = SampleID, color = IsOutlier), size = 3, vjust = -0.6) +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
    labs(
      title = paste0("cmdscale_", TaskName, " (", normalize_method, ")"),
      x = "PCoA 1", y = "PCoA 2"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  # Save PDF
  ggsave(pdf_name, plot = p, width = 12, height = 5)

  # Step 5: Update the object
  CrcBiomeScreenObject$OrginalNormalizedData <- study_data
  CrcBiomeScreenObject$NormalizedData <- study_data[!is_outlier, , drop = FALSE]
  CrcBiomeScreenObject$SampleData <- CrcBiomeScreenObject$SampleData[!is_outlier, , drop = FALSE]
  CrcBiomeScreenObject$OutlierSamples <- outliers

  return(CrcBiomeScreenObject)
}

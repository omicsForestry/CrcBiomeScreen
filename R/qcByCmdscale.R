#' Quality control using classical MDS and outlier detection
#'
#' This function performs quality control on microbiome data by applying classical
#' multidimensional scaling (MDS) on the relative abundance matrix. It identifies outlier samples
#' based on their Euclidean distance to the centroid in the first two MDS dimensions.
#'
#' Outlier samples are removed from the normalized matrix and sample metadata, and the results
#' are visualized and saved as a PDF.
#'
#' @param CrcBiomeScreenObject A \code{CrcBiomeScreenObject} containing normalized microbiome data, sample metadata, etc.
#' @param TaskName A character string used to label the output plot and PDF filename.
#' @param normalize_method A character string indicating the normalization method used (e.g., \code{"GMPR"}). Used for labeling only.
#' @param threshold_sd Numeric value indicating how many standard deviations above the mean
#'        distance should be considered an outlier (default is 1).
#'
#' @importFrom dplyr mutate across
#' @importFrom tibble tibble
#'
#'
#' @return A modified \code{CrcBiomeScreenObject} where:
#' \itemize{
#'   \item \code{$NormalizedData} contains filtered data with outliers removed.
#'   \item \code{$SampleData} is updated to exclude outlier samples.
#'   \item \code{$OutlierSamples} is a character vector of sample IDs identified as outliers.
#'   \item \code{$OrginalNormalizedData} stores the unfiltered data matrix before QC.
#' }
#'
#' @details
#' The function calculates the Euclidean distance between samples in the 2D MDS space.
#' Samples whose distance to the centroid exceeds the threshold (mean + \code{threshold_sd} * SD)
#' are considered outliers.
#'
#' A PDF plot is saved to the working directory, showing sample positions in MDS space
#' with outliers highlighted in red.
#'
#' @export
#'
#' @examples
#' # Perform QC with threshold = 1 SD
#' CrcBiomeScreenObject <- qcByCmdscale(CrcBiomeScreenObject,
#'   TaskName = "Normalize_ToyData_filtered_qc",
#'   normalize_method = "GMPR"
#' )
qcByCmdscale <- function(CrcBiomeScreenObject,
                         TaskName = NULL,
                         normalize_method = NULL,
                         threshold_sd = 1,
                         plot = TRUE) {
  # Extract normalized data and sample IDs
  study_data <- CrcBiomeScreenObject$NormalizedData
  sampleID <- rownames(CrcBiomeScreenObject$SampleData)

  # Step 1: Compute Euclidean distance matrix and apply classical MDS
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
  
  if(!is.null(outliers)){
    message(paste("Number of outliers detected:", length(outliers)))
    # Ensure IsOutlier is a factor for color mapping
    plot_df <- data.frame(
      SampleID = sampleID,
      Dim1 = x,
      Dim2 = y,
      IsOutlier = factor(is_outlier, levels = c(FALSE, TRUE))
    )
    plot_df$IsOutlier <- factor(plot_df$IsOutlier, levels = c(FALSE, TRUE))
    
    if(plot == TRUE){
    # Output file
    pdf_name <- paste0("cmdscale_", TaskName, "_", normalize_method, ".pdf")
    
    # Create plot
    p <- ggplot2::ggplot(plot_df, aes(x = Dim1, y = Dim2)) +
      geom_point(color = "grey50", size = 2) +
      geom_point(
        data = subset(plot_df, IsOutlier == TRUE),
        ggplot2::aes(x = Dim1, y = Dim2),
        color = "red", size = 2, shape = 1, stroke = 1.2
      ) +
      geom_text(aes(label = SampleID, color = IsOutlier), size = 3, vjust = -0.6) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
      labs(
        title = paste0("cmdscale_", TaskName, " (", normalize_method, ")"),
        x = "PCoA 1", y = "PCoA 2"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Save PDF
    ggsave(pdf_name, plot = p, width = 12, height = 5)
    }
  } else {
    message("No outliers detected.")
  }
  
  # Step 5: Update the object
  CrcBiomeScreenObject$OrginalNormalizedData <- study_data
  CrcBiomeScreenObject$NormalizedData <- study_data[!is_outlier, , drop = FALSE]
  CrcBiomeScreenObject$SampleData <- CrcBiomeScreenObject$SampleData[!is_outlier, , drop = FALSE]
  if(!is.null(outliers)){
    attr(outliers, "QC") <- TaskName
    attr(outliers, "OutlierSamples") <- length(outliers)
    CrcBiomeScreenObject$OutlierSamples <- outliers
  }

  return(CrcBiomeScreenObject)
}

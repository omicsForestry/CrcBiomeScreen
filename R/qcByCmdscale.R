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
#' @param plot Logical value indicating whether to generate and save the MDS plot (default is TRUE).
#' @param outdir The output directory where plots will be saved (default: tempdir()).description
#' @importFrom dplyr mutate across
#' @importFrom tibble tibble
#'
#'
#' @return A modified \code{CrcBiomeScreenObject} where:
#' \itemize{
#'   \item \code{NormalizedData} contains filtered data with outliers removed.
#'   \item \code{SampleData} is updated to exclude outlier samples.
#'   \item \code{OutlierSamples} is a character vector of sample IDs identified as outliers.
#'   \item \code{OriginalNormalizedData} stores the unfiltered data matrix before QC.
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
#' @return A A \code{CrcBiomeScreen} object. with outliers.
#' @export
#'
#' @examples
#' # Minimal toy object for QC example
#' toy_sampledata <- data.frame(
#'   sample_id = paste0("S", 1:4),
#'   study_condition = c("control", "CRC", "control", "CRC"),
#'   row.names = paste0("S", 1:4)
#' )
#'
#' # Samples in rows and features in columns
#' toy_norm <- data.frame(
#'   g1 = c(1, 2, 1, 3),
#'   g2 = c(2, 3, 1, 2),
#'   g3 = c(3, 4, 1, 1),
#'   row.names = paste0("S", 1:4)
#' )
#'
#' toy_obj <- new(
#'   "CrcBiomeScreen",
#'   AbsoluteAbundance = data.frame(),
#'   RelativeAbundance = data.frame(),
#'   TaxaData = data.frame(),
#'   SampleData = toy_sampledata,
#'   NormalizedData = toy_norm,
#'   ModelData = list()
#' )
#'
#' qc_obj <- qcByCmdscale(
#'   toy_obj,
#'   TaskName = "ToyQC",
#'   normalize_method = "GMPR",
#'   threshold = 1
#' )
#
qcByCmdscale <- function(CrcBiomeScreenObject,
                         TaskName = NULL,
                         outdir = tempdir(),
                         normalize_method = NULL,
                         threshold_sd = 1,
                         plot = FALSE) {
  # Extract normalized data and sample IDs
  study_data <- CrcBiomeScreenObject@NormalizedData
  sampleID <- rownames(CrcBiomeScreenObject@SampleData)

  # Safety check: samples should be arranged in rows
  if (nrow(study_data) != nrow(CrcBiomeScreenObject@SampleData)) {
    stop(
      "The number of rows in NormalizedData does not match the number of rows ",
      "in SampleData. qcByCmdscale() expects samples in rows and features in columns."
    )
  }

  if (length(sampleID) != nrow(study_data)) {
    stop(
      "The length of sampleID does not match the number of rows in NormalizedData."
    )
  }

  # Step 1: Compute Euclidean distance matrix and apply classical MDS
  dist_data <- stats::dist(study_data)
  mds_coords <- stats::cmdscale(dist_data, k = 2)
  rownames(mds_coords) <- sampleID

  x <- mds_coords[, 1]
  y <- -mds_coords[, 2] # Flip axis 2 for better visualization

  # Step 2: Compute Euclidean distance of each sample to the centroid
  center <- c(mean(x), mean(y))
  distances <- sqrt((x - center[1])^2 + (y - center[2])^2)

  # Step 3: Identify outliers based on mean + n*SD rule
  threshold <- mean(distances) + threshold_sd * stats::sd(distances)
  is_outlier <- distances > threshold
  outliers <- sampleID[is_outlier]

  # Safety checks before removing samples
  if (length(is_outlier) != nrow(study_data)) {
    stop(
      "Length of outlier indicator does not match the number of rows in ",
      "NormalizedData. Please check whether samples are arranged in rows."
    )
  }

  if (length(is_outlier) != nrow(CrcBiomeScreenObject@SampleData)) {
    stop(
      "Length of outlier indicator does not match the number of rows in SampleData."
    )
  }

  # Print outlier information
  message("Outlier sample IDs detected:")

  if (length(outliers) > 0) {
    message(paste(outliers, collapse = "\n"))
    message(sprintf("Number of outliers detected: %d", length(outliers)))
  } else {
    message("No outliers detected.")
  }

  # Prepare plot data
  plot_df <- data.frame(
    SampleID = sampleID,
    Dim1 = x,
    Dim2 = y,
    IsOutlier = factor(is_outlier, levels = c(FALSE, TRUE))
  )

  # Circle data for threshold boundary
  circle_df <- data.frame(
    x = center[1] + threshold * cos(seq(0, 2 * pi, length.out = 200)),
    y = center[2] + threshold * sin(seq(0, 2 * pi, length.out = 200))
  )

  # Plot regardless of whether outliers are detected
  if (isTRUE(plot)) {
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

    pdf_name <- paste0("cmdscale_", TaskName, "_", normalize_method, ".pdf")

    plot_subtitle <- paste0(
      "Threshold = mean + ", threshold_sd, " SD; ",
      "Outliers detected: ", length(outliers)
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2)) +
      ggplot2::geom_point(
        ggplot2::aes(color = IsOutlier),
        size = 2
      ) +
      ggplot2::geom_path(
        data = circle_df,
        ggplot2::aes(x = x, y = y),
        color = "blue",
        linetype = "dashed",
        linewidth = 0.8,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_point(
        data = data.frame(x = center[1], y = center[2]),
        ggplot2::aes(x = x, y = y),
        color = "blue",
        size = 2,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = c(`FALSE` = "grey50", `TRUE` = "red"),
        drop = FALSE
      ) +
      ggplot2::labs(
        title = paste0("cmdscale_", TaskName, " (", normalize_method, ")"),
        subtitle = plot_subtitle,
        x = "MDS 1",
        y = "MDS 2",
        color = "Outlier"
      ) +
      ggplot2::theme_minimal()

    # Label only outliers if any are detected
    if (length(outliers) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = subset(plot_df, IsOutlier == TRUE),
          ggplot2::aes(x = Dim1, y = Dim2),
          color = "red",
          size = 2,
          shape = 1,
          stroke = 1.2
        ) +
        ggplot2::geom_text(
          data = subset(plot_df, IsOutlier == TRUE),
          ggplot2::aes(label = SampleID),
          color = "red",
          size = 3,
          vjust = -0.6
        )
    }

    ggplot2::ggsave(
      filename = file.path(outdir, pdf_name),
      plot = p,
      width = 12,
      height = 5
    )
  }

  # Store original normalized data
  CrcBiomeScreenObject@OriginalNormalizedData <- study_data

  # Remove outliers only if detected
  if (length(outliers) > 0) {
    CrcBiomeScreenObject@NormalizedData <- study_data[!is_outlier, , drop = FALSE]
    CrcBiomeScreenObject@SampleData <- CrcBiomeScreenObject@SampleData[!is_outlier, , drop = FALSE]

    attr(outliers, "QC") <- TaskName
    attr(outliers, "OutlierSamples") <- length(outliers)

    CrcBiomeScreenObject@OutlierSamples <- outliers
  } else {
    CrcBiomeScreenObject@NormalizedData <- study_data
    CrcBiomeScreenObject@OutlierSamples <- character(0)
  }

  return(CrcBiomeScreenObject)
}

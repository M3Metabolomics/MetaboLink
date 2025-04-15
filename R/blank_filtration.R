#' Blank Filtration: filter data based on blank and QC sample sigal strength 
#' 
#' @param data A data.frame or matrix with signal intensities.
#' @param sequence A matrix or data.frame indicating sample types and metadata.
#' @param signalStrength A numeric multiplier to threshold blank signals.
#' @param keepIs A logical indicating whether to keep internal standards (IS).
#' @return A filtered version of the input data.
#' 
blankFiltration <- function(data, sequence, signalStrength, keepIs) {
  # replace NAs in blank columns with 0
  data[sequence[, 1] %in% "Blank"][is.na(data[sequence[, 1] %in% "Blank"])] <- 0

  # identify features when mean blank * signal strength is less than the mean QC signal
  keepFeature <- apply(data[sequence[, 1] %in% "Blank"], 1, mean) * signalStrength < 
                    apply(data[sequence[, 1] %in% "QC"], 1, mean, na.rm = TRUE)

  if (keepIs) {
    # keep features that are internal standards (IS)
    isFeature <- grepl("\\(IS\\)", toupper(data[sequence[, 1] %in% "Name"][, 1]))
    data <- data[bf | isFeature, ]
  } else {
    data <- data[bf, ]
  }

  return(data)
}
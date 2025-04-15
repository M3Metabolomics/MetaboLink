
#' Impute Missing Values by Group
#'
#' @details The function combines the input data with the grouping variable and
#'   transposes the data for easier manipulation. It then groups the data by
#'   the `seq$group` variable and replaces missing values in numeric columns
#'   based on the specified `method`.
#'
impute_missing_by_group <- function(data, seq, method = c("Median", "Min/X"), minx = NULL) {
  method <- match.arg(method)

  datm <- data.frame(group = seq$group, t(data))

  datm <- datm %>%
    dplyr::group_by(group) %>%
    dplyr::mutate_if(
      is.numeric,
      function(x) {
        if (method == "Median") {
          ifelse(is.na(x), median(x, na.rm = TRUE), x)
        } else if (method == "Min/X") {
          if (is.null(minx)) {
            stop("Argument 'minx' must be provided when method is 'minx'.")
          }
          min_val <- min(x, na.rm = TRUE)
          ifelse(is.na(x), min_val / minx, x)
        } else {
          # This should not happen due to match.arg
          stop("Invalid imputation method.")
        }
      }
    )
  
  # transpose back and remove group column
  datm <- data.frame(t(datm[, -1])) 

  # restore column and row names
  colnames(datm) <- colnames(data)
  row.names(datm) <- row.names(data)

  datm[datm == "Inf"] <- NA

  return(datm)
}


#' Imputation Function
#'
#' @description Performs missing value imputation on a given dataset based on the specified method.
#'
imputation <- function(data, seq, method, minx = 1, process_only_qc, remaining_strategy) {
  data_copy <- data

  if (process_only_qc) {
    qcData <- data[seq[, 1] %in% "QC"]
    qcSequence <- seq[seq[, 1] %in% "QC", ]
  } else {
    qcData <- data[seq[, 1] %in% c("Sample", "QC")]
    qcSequence <- seq[seq[, 1] %in% c("Sample", "QC"), ]
  }

  qcData[qcData == 0] <- NA
  qcSequence[qcSequence[, 1] %in% c("QC"), ]$group <- "QC"
  qcSequence$group[is.na(qcSequence$group)] <- "Sample"

  if (method == "KNN") {
    qcData <- as.matrix(qcData)
    knndat <- impute.knn(qcData, k = 10, rowmax = .99, colmax = .99, maxp = 15000) # TODO should k be set differently, e.g., smallest class -1?
    impsqdat <- as.data.frame(knndat$data)

  } else {
    impsqdat <- impute_missing_by_group(qcData, qcSequence, method = method, minx = minx)
  }

  # if there's still missing values after imputation, apply the remaining strategy
  if (sum(is.na(impsqdat)) > 0) {
    if (remaining_strategy == "Min/X") {
      for (i in 1:nrow(impsqdat)) {
        impsqdat[i, is.na(impsqdat[i, ])] <- min(impsqdat[i, ], na.rm = T) / minx
      }
    } else if (remaining_strategy == "zero") {
      for (i in 1:nrow(impsqdat)) {
        impsqdat[i, is.na(impsqdat[i, ])] <- 0
      }
    } else if (remaining_strategy == "Median") {
      for (i in 1:nrow(impsqdat)) {
        impsqdat[i, is.na(impsqdat[i, ])] <- median(as.numeric(impsqdat[i, ]), na.rm = T)
      }
    }
  }


  # replace imputed values in the original data
  if (process_only_qc) {
    data_copy[seq[, 1] %in% "QC"] <- impsqdat
  } else {
    data_copy[seq[, 1] %in% c("Sample", "QC")] <- impsqdat
  }
  
  return(data_copy)
}
#' Normalization Using Internal Standards
#'
#' @details
#' The function normalizes the data by dividing each value by the corresponding value of the nearest internal standard. If the method is "Same lipid structure", the function ensures that the internal standard used for normalization has the same lipid structure as the target compound.
#'
normalizationIS <- function(data, sequence, is, method, qc) {

  sendSweetAlert(title = "Normalizing...", text = "This modal will close automatically when done. Please wait...", type = "info")

  rt_index <- which(sequence[, 1] == "RT")
  isname <- is
  is <- as.numeric(gsub(" .*$", "", is))

  sel <- if (qc) c("Sample", "QC") else "Sample"
  sdat <- data[sequence[, 1] %in% sel]

  sdat[sdat == 0] <- NA
  is <- is[complete.cases(sdat[is, ])]

  if(length(is) == 0) {
    closeSweetAlert()
    return(NULL)
  }

  near <- sapply(data[, rt_index], function(y) {
    which.min(abs(data[is, rt_index] - y))
  })

  # if the method is "Same lipid structure", find the nearest internal standard with the same lipid structure
  if (method == "Same lipid structure") {
    name <- data[sequence[, 1] %in% "Name"]
    istype <- gsub(" .*$", "", name[is, ])
    near <- sapply(seq(name[, 1]), function(x) {
      if (gsub(" .*$", "", name[x, 1]) %in% istype) {
        which(istype %in% gsub(" .*$", "", name[x, 1]))
      } else {
        near[x]
      }
    })
  }

  # normalize the data
  sdat <- sapply(seq(ncol(sdat)), function(j) {
    sapply(seq(nrow(sdat)), function(i) {
      sdat[i, j] <- sdat[i, j] / sdat[is, j][near[i]]
    })
  })

  # save the name of the internal standard used for normalization
  isnorm <- sapply(seq(nrow(sdat)), function(x) {
    isname[near[x]]
  })

  data[sequence[, 1] %in% sel] <- sdat
  data <- cbind(data, data.frame(isnorm = isnorm))
  closeSweetAlert()
  return(data)
}

findInternalStandards <- function(data) {
  # Check if the first column of the data contains internal standards
  isIndex <- grepl("\\(IS\\)", toupper(data[, 1]))

  if (sum(isIndex) > 0) {
    featureName <- as.vector(data[isIndex, 1])
    # combine indices and names
    internalStandards <- paste(which(isIndex), " - ", featureName)

    return(internalStandards)
  }
  # no IS found
  return(character(0))
}

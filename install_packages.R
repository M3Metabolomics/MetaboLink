# Install required packages for MetaboLink

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# CRAN packages
cran_packages <- c(
  "shiny",
  "shinydashboard", 
  "shinyBS",
  "shinyjs",
  "shinyWidgets",
  "spsComps",
  "DT",
  "dplyr",
  "plotly",
  "ggplot2",
  "ggrepel",
  "gridExtra",
  "randomForest",
  "writexl",
  "igraph",
  "stringi",
  "shinycssloaders",
  "jsonlite",
  "shinyalert",
  "shinybusy",
  "gtools",
  "caret",
  "ggbeeswarm",
  "broom",
  "webchem",
  "PubChemR",
  "data.table",
  "car",
  "stringr",
  "circlize",
  "grid",
  "stats",
  "colourpicker",
  "scales",
  "ggraph",
  "tidygraph",
  "tidyr",
  "networkD3",
  "network",
  "sna",
  "visNetwork",
  "threejs",
  "ndtv",
  "ggnetwork",
  "devtools",
  "factoextra",
  "ggdendro",
  "dbscan",
  "PubChemR"
)

# Bioconductor packages
bioc_packages <- c(
  "impute",
  "ComplexHeatmap",
  "InteractiveComplexHeatmap",
  "clusterProfiler",
  "KEGGREST",
  "rgoslin",
  "lipidomeR"
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
install.packages(cran_packages, dependencies = TRUE)

# Install Bioconductor packages
cat("Installing Bioconductor packages...\n")
BiocManager::install(bioc_packages, dependencies = TRUE)

# Install PubChemR from GitHub if not already installed
if (!requireNamespace("PubChemR", quietly = TRUE)) {
  devtools::install_github("selcukorkmaz/PubChemR")
}
# Verify installation
cat("Verifying package installation...\n")
all_packages <- c(cran_packages, bioc_packages, "BiocManager")

for (pkg in all_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "installed successfully\n")
  } else {
    cat("✗", pkg, "installation failed\n")
  }
}

cat("Package installation complete!\n")
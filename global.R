# Base
library(grid)
library(stats)

# Shiny 
library(shiny)
library(shinyalert)
library(shinyBS)
library(shinybusy) # for spinners, progress bar...
library(shinycssloaders)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(spsComps) # custom  UI and server components
library(colourpicker)

# Data science
library(tidyverse)
library(data.table) # extension of data.frame
library(DT)
library(impute)
library(plotly)

library(car)
library(caret)
library(clusterProfiler)
library(randomForest)

library(circlize)
library(ggbeeswarm) # categorical scatter (violin point) plots
library(ggnetwork) # ggplot2 extension for network visualization
library(ggrepel) #base improved text labels in plots for ggplot2
library(ggraph)
library(igraph)
library(tidygraph)

library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library("lipidomeR")

library(ndtv)
library(network)
library(networkD3)
library(sna)
library(threejs)
library(visNetwork)

# Utils
library(BiocManager)
library(broom)
library(gridExtra)
library(jsonlite)
library(rgoslin) # lipid name parser
library(scales)
library(stringi)
library(writexl)

# APIs
library(KEGGREST)
library(PubChemR)
library(webchem)


# library(gtools) # R programming tools


source("functions.R")

# Source files in R folder
rFiles <- list.files("./R", pattern = "\\.R$", full.names = TRUE)
for (file in rFiles) {
  source(file)
}
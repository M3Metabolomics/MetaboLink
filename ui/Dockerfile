FROM rocker/shiny
LABEL maintainer="Ana Mendes <anamendesml@outlook.com>"
LABEL description="Docker image of MetaboLink"

RUN rm -rf /srv/shiny-server/*
COPY . /srv/shiny-server/
WORKDIR /srv/shiny-server/
RUN apt update; apt install -y libglpk-dev

# RUN apt-get -o Acquire::AllowInsecureRepositories=true -o Acquire::AllowDowngradeToInsecureRepositories=true update && apt-get install -y libssl-dev liblzma-dev libbz2-dev libicu-dev && apt-get clean

#### CRAN ####
# Shiny
RUN R -e "install.packages(c('shiny', 'shinyalert', 'shinyBS', 'shinybusy', 'shinycssloaders', 'shinydashboard',  'shinyjs',  'shinyWidgets', 'spsComps', 'colourpicker'), repos='http://cran.us.r-project.org')"

# Data science, visualization
RUN R -e "install.packages(c('tidyverse', 'data.table', 'DT', 'plotly', 'car', 'caret', 'randomForest'), repos='http://cran.us.r-project.org')"

RUN R -e "install.packages(c('circlize', 'ggbeeswarm', 'ggnetwork', 'ggrepel', 'ggraph', 'igraph', 'tidygraph'), repos='http://cran.us.r-project.org')"

RUN R -e "install.packages(c('ndtv', 'network', 'networkD3', 'sna', 'threejs', 'visNetwork', 'lipidomeR'), repos='http://cran.us.r-project.org')"

# Utils
RUN R -e "install.packages(c('BiocManager', 'broom', 'gridExtra', 'jsonlite', 'scales', 'stringi', 'writexl'), repos='http://cran.us.r-project.org')"

RUN R -e "install.packages(c('PubChemR', 'webchem'), repos='http://cran.us.r-project.org')"


#### Bioconductor ####
RUN R -e "library(BiocManager); BiocManager::install(c('impute', 'rgoslin', 'KEGGREST', 'limma', 'qvalue', 'ComplexHeatmap', 'InteractiveComplexHeatmap', 'clusterProfiler'), ask=F)"

# tibble, ggdentro, dbscan, factoextra, htmltools
# 'enrichplot'

# ggiraph, GO.db
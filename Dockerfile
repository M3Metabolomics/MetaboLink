FROM rocker/shiny
LABEL maintainer="Ana Mendes <anamendesml@outlook.com>"
LABEL description="Docker image of MetaboLink"

RUN rm -rf /srv/shiny-server/*
COPY . /srv/shiny-server/
WORKDIR /srv/shiny-server/

RUN apt-get update && \
    apt-get install -y --no-install-recommends libglpk-dev libmagick++-dev imagemagick && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript /srv/shiny-server/install_packages.R

RUN R -e "library(devtools); devtools::install_github('selcukorkmaz/PubChemR')"

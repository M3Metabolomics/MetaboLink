FROM rocker/geospatial:4.5.1
RUN apt-get update -y && apt-get install -y  libcurl4-openssl-dev libpng-dev libssl-dev perl make libicu-dev pandoc cmake libfontconfig1-dev libfreetype6-dev libfribidi-dev libharfbuzz-dev libxml2-dev zlib1g-dev libmagick++-dev gsfonts librsvg2-dev texlive libx11-dev libglpk-dev libjpeg-dev libtiff-dev libwebp-dev git && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(renv.config.pak.enabled = FALSE, repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_version("renv", version = "1.0.3")'
COPY renv.lock renv.lock
RUN --mount=type=cache,id=renv-cache,target=/root/.cache/R/renv R -e 'renv::restore()'
WORKDIR /srv/shiny-server/
COPY . /srv/shiny-server/
EXPOSE 3838
CMD R -e 'shiny::runApp("/srv/shiny-server",host="0.0.0.0",port=3838)'

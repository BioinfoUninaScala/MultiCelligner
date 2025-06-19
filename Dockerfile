FROM rocker/shiny:latest

LABEL maintainer="Giovanni Scala giovanni.scala@unina.it"

# system libraries of general use
RUN apt-get update && apt-get install -y libssl-dev libxml2-dev libfontconfig1-dev libcurl4-openssl-dev libharfbuzz-dev libfribidi-dev libpq-dev libglpk-dev git

# R packages

RUN R -e "install.packages(c('remotes','devtools', 'fontawesome', 'htmltools','ggplot2','SNFtool','reactable','shiny', 'plotly', 'dplyr','stringr', 'reshape2','crosstalk', 'magrittr', 'shinyjs', 'shinycssloaders'))" 
RUN R -e "remotes::install_github('BioinfoUninaScala/MultiCelligner', build_vignettes=FALSE,dependencies=TRUE, type='source',upgrade = 'never')"

COPY Rprofile.site /usr/local/lib/R/etc/

EXPOSE 3838
CMD ["R", "-e", "shiny::runApp(MultiCelligner::MultiCellignerShiny(), host = '0.0.0.0', port = 3838)"]
FROM rocker/shiny:latest

LABEL maintainer="Giovanni Scala giovanni.scala@unina.it"

# system libraries of general use
RUN apt-get update && apt-get install -y libssl-dev libxml2-dev libfontconfig1-dev libcurl4-openssl-dev libharfbuzz-dev libfribidi-dev libpq-dev libglpk-dev git



# R packages

RUN R -e "install.packages(c('remotes','devtools', 'BiocManager','shinyjs'))" && R -e "BiocManager::install(version = '3.19',ask = FALSE); BiocManager::install(c('S4Arrays', 'XVector', 'SparseArray','rhdf5','S4Vectors','assorthead','HDF5Array','BiocNeighbors', 'BiocParallel', 'maftools', 'MOFA2'))"

# RUN R -e "remotes::install_github('BioinfoUninaScala/MultiCelligner', build_vignettes=FALSE,dependencies=TRUE, type='source',upgrade = "never")"
# RUN R -e "remotes::install_github('daattali/shinycssloaders')"

COPY MultiCelligner /MultiCelligner
RUN R -e "devtools::install('/MultiCelligner/',dependencies = TRUE, upgrade = 'never')"

COPY Rprofile.site /usr/local/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "MultiCelligner::MultiCellignerShiny()"]
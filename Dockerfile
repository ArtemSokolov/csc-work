FROM rocker/tidyverse:latest

RUN R -e 'install.packages("languageserver", "egg")'
RUN R -e 'install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))'
RUN R -e 'devtools::install_github("ArtemSokolov/synExtra")'
RUN R -e 'devtools::install_github("labsyspharm/naivestates")'

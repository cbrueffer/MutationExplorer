FROM rocker/shiny:3.5.3
LABEL maintainer="Christian Brueffer <christian.brueffer@med.lu.se>"


###############################################
#
# General
#
###############################################

EXPOSE 3838/tcp


###############################################
#
# Package installation
#
###############################################

# install and configure the nginx shiny proxy
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
                       libssl-dev \
                       libbz2-dev \
                       liblzma-dev \
                       libxml2-dev

RUN R -e "install.packages(c( \
    'cowplot', \
    'dbplyr', \
    'dplyr', \
    'devtools', \
    'DT', \
    'DBI', \
    'ggrepel', \
    'httr', \
    'jsonlite', \
    'reshape2', \
    'RColorBrewer', \
    'RSQLite', \
    'sessioninfo', \
    'shinycssloaders', \
    'shinyhelper', \
    'shinyjs', \
    'shinyWidgets', \
    'stringr', \
    'survminer', \
    'yaml', \
    'BiocManager' \
    ), clean = TRUE)"
RUN R -e "library(devtools); install_github('cbrueffer/survminer', ref='arrange_flexibility')"
RUN R -e "BiocManager::install(c( \
    'reactome.db', \
    'GenVisR' \
    ), clean = TRUE, version = '3.8', ask = FALSE, update = TRUE)"

###############################################
#
# Application installation
#
###############################################

ENV APP_LOCATION /srv/shiny-server/MutationExplorer

COPY *.R ${APP_LOCATION}/
COPY config.yaml ${APP_LOCATION}/
COPY google-analytics.js ${APP_LOCATION}/
COPY about.md ${APP_LOCATION}/
COPY R ${APP_LOCATION}/R
COPY helpfiles ${APP_LOCATION}/helpfiles
COPY index.html /srv/shiny-server/index.html

CMD /usr/bin/shiny-server.sh

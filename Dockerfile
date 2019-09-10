FROM rocker/shiny:latest
LABEL maintainer="Christian Brueffer <christian.brueffer@med.lu.se>"


###############################################
#
# General
#
###############################################

EXPOSE 80/tcp
# shiny-server not exposed; accessible through the nginx proxy
#EXPOSE 3838/tcp


###############################################
#
# Package installation
#
###############################################

# install and configure the nginx shiny proxy
RUN apt-get update && \
    apt-get install -y nginx \
                       libssl-dev \
                       libxml2-dev
COPY nginx-shiny-proxy.conf /etc/nginx/nginx.conf

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
    'magrittr', \
    'reshape2', \
    'RColorBrewer', \
    'RSQLite', \
    'shinycssloaders', \
    'shinyhelper', \
    'shinyjs', \
    'shinyWidgets', \
    'stringr', \
    'survminer', \
    'yaml', \
    'BiocManager' \
    ))"
RUN R -e "library(devtools); install_github('cbrueffer/survminer', ref='arrange_flexibility')"
RUN R -e "BiocManager::install(c( \
    'reactome.db', \
    'GenVisR' \
    ), version = '3.8', ask = FALSE, update = TRUE)"

###############################################
#
# Application installation
#
###############################################

ENV APP_LOCATION /srv/shiny-server/MutationExplorer

COPY *.R ${APP_LOCATION}/
COPY config.yaml ${APP_LOCATION}/
COPY htpasswd-sme.txt ${APP_LOCATION}/
COPY about.md ${APP_LOCATION}/
COPY R ${APP_LOCATION}/R
COPY data ${APP_LOCATION}/data
COPY helpfiles ${APP_LOCATION}/helpfiles

CMD /etc/init.d/nginx start && /usr/bin/shiny-server.sh

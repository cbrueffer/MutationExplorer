FROM rocker/shiny:latest
LABEL maintainer="Christian Brueffer <christian.brueffer@med.lu.se>"


###############################################
#
# General
#
###############################################

EXPOSE 3838/tcp
EXPOSE 80/tcp


###############################################
#
# Package installation
#
###############################################

RUN apt-get update && apt-get install -y nginx
COPY mutationexplorer.conf /etc/nginx/conf.d/mutationexplorer.conf

RUN R -e "install.packages(c('dplyr', 'DT', 'magrittr', 'shinycssloaders', 'shinyhelper', 'shinyjs', 'shinyWidgets', 'survminer', 'BiocManager'))"
RUN R -e "BiocManager::install('reactome.db', version = '3.8', ask = FALSE, update = TRUE)"


###############################################
#
# Application installation
#
###############################################

ENV APP_LOCATION /srv/shiny-server/ShinyMutationExplorer

COPY *.R ${APP_LOCATION}/
COPY R ${APP_LOCATION}/R
COPY data ${APP_LOCATION}/data
COPY helpfiles ${APP_LOCATION}/helpfiles

CMD ["/usr/bin/shiny-server.sh"]

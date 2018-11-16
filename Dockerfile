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

RUN R -e "install.packages(c('dplyr', 'DT', 'knitr', 'magrittr', 'shinycssloaders', 'shinyhelper', 'shinyjs', 'shinyWidgets', 'survminer', 'BiocManager'))"
RUN R -e "BiocManager::install('reactome.db', version = '3.8', ask = FALSE, update = TRUE)"

# install and configure the nginx shiny proxy
RUN apt-get update && \
    apt-get install -y nginx
COPY nginx-shiny-proxy.conf /etc/nginx/nginx.conf

###############################################
#
# Application installation
#
###############################################

ENV APP_LOCATION /srv/shiny-server/ShinyMutationExplorer

COPY *.R ${APP_LOCATION}/
COPY htpasswd-sme.txt ${APP_LOCATION}/
COPY R ${APP_LOCATION}/R
COPY data ${APP_LOCATION}/data
COPY helpfiles ${APP_LOCATION}/helpfiles

CMD /etc/init.d/nginx start && /usr/bin/shiny-server.sh

#FROM rstudio/r-base:3.5-xenial
FROM jupyter/r-notebook
USER root

WORKDIR app

RUN sudo apt-get update && \
    sudo apt-get upgrade -y && \
    sudo apt-get install -y git

RUN git clone https://github.com/erenasena/bm-simulations.git

RUN R -e "install.packages(c('fExpressCertificates', 'somebm', 'yuima'), repo='http://cran.rstudio.com/')"

EXPOSE 8888

CMD jupyter notebook --allow-root

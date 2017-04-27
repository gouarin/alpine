FROM frolvlad/alpine-glibc

MAINTAINER Loic Gouarin "loic.gouarin@gmail.com"

# Configure environment
ENV CONDA_DIR=/opt/conda CONDA_VER=4.3.11
ENV PATH=$CONDA_DIR/bin:$PATH SHELL=/bin/bash LANG=C.UTF-8

# Install useful packages 
RUN echo "@testing http://nl.alpinelinux.org/alpine/edge/testing" >> /etc/apk/repositories &&\
    apk --update add \
    bash \
    curl \
    git \
    ca-certificates \
    tini@testing \
    libice \
    libsm \
    libstdc++ \
    freeglut &&\
    rm -rf /var/cache/apk/*

RUN curl https://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-x86_64.sh  -o mconda.sh 

RUN /bin/bash mconda.sh -f -b -p $CONDA_DIR && \
    rm mconda.sh
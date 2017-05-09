FROM debian:sid-slim

MAINTAINER Loic Gouarin "loic.gouarin@gmail.com"

# Configure environment
ENV USER=precis

# Install useful packages 
RUN apt-get update && \
    apt-get install -y curl bzip2 build-essential python libopenblas-dev cmake openssh-server && \
    apt-get clean

# get and install miniconda
RUN curl https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.0.tar.gz -o openmpi-2.1.0.tar.gz && \
    tar xzf openmpi-2.1.0.tar.gz && \
    cd openmpi-2.1.0 && \
    ./configure \
            --disable-dependency-tracking \
            --enable-mpi-cxx && \
    make install && \
    cd / && rm -r openmpi-2.1.0.tar.gz openmpi-2.1.0

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

RUN curl http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz -o petsc.tar.gz && \
    tar xzf petsc.tar.gz && \
    cd petsc-3.7.6 && \
    CC=mpicc CXX=mpicxx python ./configure \
        --with-fc=0 \
        --with-debugging=0 \
        --COPTFLAGS=-O3 \
        --CXXOPTFLAGS=-O3 \
        --with-hwloc=0 \
        --with-mpi=1 \
        --with-pthread=1 \
        --with-ssl=0 \
        --with-suitesparse=0 \
        --with-x=0 && \
    make PETSC_DIR=/petsc-3.7.6 PETSC_ARCH=arch-linux2-c-opt all && \
    rm /petsc.tar.gz

ENV PETSC_DIR=/petsc-3.7.6
ENV PETSC_ARCH=arch-linux2-c-opt

RUN useradd -m $USER

WORKDIR /home/$USER

COPY Poisson /home/$USER/Poisson
RUN chown -R $USER:$USER /home/$USER/Poisson

USER $USER


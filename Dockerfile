FROM ubuntu:18.04

RUN apt-get update \
    && apt-get install -y git dpkg-dev make g++ gcc binutils libx11-dev libxpm-dev \
    libxft-dev libxext-dev gfortran libssl-dev libpcre3-dev \
    xlibmesa-glu-dev libglew1.5-dev libftgl-dev \
    libmysqlclient-dev libfftw3-dev \
    graphviz-dev libavahi-compat-libdnssd-dev \
    libldap2-dev python-dev libxml2-dev libkrb5-dev \
    libgsl0-dev libqt4-dev python-pip wget

WORKDIR /root

RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz && \
    tar -xzvf jellyfish-2.2.10.tar.gz && \
    cd jellyfish-2.2.10 && \
    ./configure --prefix=/root/jellyfish-2.2.10 && \
    make && \
    make install && \
    rm ../jellyfish-2.2.10.tar.gz

ENV PATH "/root/jellyfish-2.2.10/bin:${PATH}"

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    rm ncbi-blast-2.7.1+-x64-linux.tar.gz

ENV PATH "/root/ncbi-blast-2.7.1+/bin:${PATH}"

COPY . /root/MHcut

RUN pip2 install -e /root/MHcut


FROM ubuntu:latest

MAINTAINER nobody@a.de
LABEL description="github actions test"

ENV ROOT_ARCHIVE="root_v6.18.04.Linux-ubuntu18-x86_64-gcc7.4.tar.gz"
ENV ROOTSYS=/opt/root

RUN apt update && apt install -y \
    build-essential bash gfortran cmake \
    libboost-all-dev wget libblas3 python-pip \
    libxml2-utils \
    doxygen libxslt1-dev xsltproc && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt && cd /opt && \
    wget https://root.cern/download/${ROOT_ARCHIVE} && \
    tar xvf ${ROOT_ARCHIVE} && \
    rm -f ${ROOT_ARCHIVE}

ADD run_action.sh /run_action.sh
ENTRYPOINT ["/bin/bash", "-x"]
CMD ["/run_action.sh"]

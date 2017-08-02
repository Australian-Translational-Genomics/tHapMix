FROM ubuntu:16.04

# install build and config tools
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    pkg-config \
    liblzma-dev \
    libbz2-dev \
    libatlas-base-dev \
    cmake \
    ncurses-dev \
    zlib1g-dev \
    software-properties-common

# install python 2.7
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    openssl \
    python2.7 \
    python2.7-dev \
    python-pip \
    python-software-properties

# install tHapMix non-pip dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    bedtools \
    samtools \
    gfortran

# copy git repository into the image
RUN mkdir -p /opt/thapmix
COPY . /opt/thapmix/
WORKDIR /opt/thapmix/

# install pip dependencies
RUN pip install --upgrade setuptools wheel
RUN pip install -r ./requirements.txt

RUN DEBIAN_FRONTEND=noninteractive apt-get clean

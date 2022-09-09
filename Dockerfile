FROM ubuntu

MAINTAINER Honzo <hz_juanito18@outlook.com>
LABEL description="Honzo's pipeline for in-silico PrimerDesign"

ARG DEBIAN_FRONTEND=noninteractive

# Get ubuntu pre-requisites
RUN apt-get update -y
RUN apt-get install -y \
        build-essential \
        python3 \
        pip \
        python-is-python3 \
        muscle \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

COPY /opt /opt
WORKDIR /opt
RUN pip install -r requirements.txt

# install argtable (required for clustalw)
WORKDIR /opt
RUN tar -xvzf argtable2-13.tar.gz && \
    cd argtable2-13 && \
    ./configure && \
    make && make install && make clean

# install clustal
WORKDIR /opt
RUN tar -xvzf clustal-omega-1.2.4.tar.gz && \
    cd clustal-omega-1.2.4 && \
    ./configure && \
    make && make install

RUN rm *.bz2 *.zip *.gz

COPY /database /database

COPY /app /app
ENV PATH="${PATH}:/app"

# Use Ubuntu 18.04 as the base image
FROM ubuntu:18.04

# Maintainer information
MAINTAINER Rick Wertenbroek <luoxiaolong2021@email.szu.edu.cn>

# Set non-interactive frontend to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND noninteractive

# Install required packages and clean up to keep the image lean
RUN apt-get update && apt-get -y upgrade && apt-get install -y \
    g++ make git zlib1g-dev libbz2-dev liblzma-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Create project directory and set it as the working directory
WORKDIR /home/Project_gsc

# Clone the GSC project and build it
RUN git clone https://github.com/luo-xiaolong/GSC.git && \
    cd GSC && \
    make clean && \
    make

# Set the working directory for the container
WORKDIR /home/Project_gsc/GSC



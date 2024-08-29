# Use Ubuntu as the base image
FROM ubuntu:24.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    wget \
    gfortran \
    python3-dev \
    python3-pip \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Set Python 3 as the default
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1


RUN mkdir -p ~/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
RUN bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
RUN rm -rf ~/miniconda3/miniconda.sh

RUN ~/miniconda3/bin/conda init bash
RUN ~/miniconda3/bin/conda init zsh

RUN ~/miniconda3/bin/conda create -n "venv" python=3.8
RUN ~/miniconda3/envs/venv/bin/python3 -m pip install numpy
RUN ~/miniconda3/envs/venv/bin/python3 -m pip install kamodo
RUN ~/miniconda3/envs/venv/bin/python3 -m pip install notebook



# Copy your project files into the container
WORKDIR /kamodo-wrapper
COPY ../kamodo_ccmc/kamodo-wrapper/ /kamodo-wrapper

# Build your project here, or run commands to set up your environment

RUN echo "export CPLUS_INCLUDE_PATH=/kamodo-wrapper/:${CPLUS_INCLUDE_PATH}" >> /root/.bashrc
RUN echo "export CMAKE_PREFIX_PATH=/kamodo-wrapper/:${CMAKE_PREFIX_PATH}" >> /root/.bashrc



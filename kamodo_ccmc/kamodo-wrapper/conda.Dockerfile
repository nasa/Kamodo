
# Use Miniforge3 as the base image
FROM condaforge/miniforge3

RUN conda install python=3.8


# Avoid interactive dialogues by apt-get and others
ENV DEBIAN_FRONTEND=noninteractive

# Update and install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    python-dev \
    && rm -rf /var/lib/apt/lists/*


# Update and install system dependencies and GFortran
RUN conda install gfortran cmake g++


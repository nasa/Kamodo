# Borrowed from https://github.com/xbonnin
# https://github.com/spacepy/spacepy/issues/530#issuecomment-859653900

FROM ccmc/kamodo_ccmc

RUN apt-get update
RUN apt-get install -y gfortran make
RUN conda install numpy scipy matplotlib networkx h5py sphinx numpydoc astropy

# IRBEM compiled by spacepy
# WORKDIR /
# RUN git clone https://github.com/asherp/IRBEM.git

# # install irbem
# WORKDIR /IRBEM
# RUN make all.help
# RUN make OS=linux64 ENV=gfortran64 all
# RUN make OS=linux64 ENV=gfortran64 install


# install spacepy
WORKDIR /
RUN git clone https://github.com/asherp/spacepy.git
WORKDIR /spacepy

# install netcdf4
RUN conda install netCDF4

RUN python setup.py install

WORKDIR /kamodo/kamodo/flythrough

# RUN python SatelliteFlythrough.py 
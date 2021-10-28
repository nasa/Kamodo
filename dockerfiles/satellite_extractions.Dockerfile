FROM ccmc/kamodo_ccmc


RUN pip install 'git+https://github.com/hapi-server/client-python' --upgrade

RUN conda install -c astropy -c defaults \
  scipy h5py beautifulsoup4 html5lib bleach pyyaml pandas sortedcontainers \
  pytz matplotlib setuptools mpmath bottleneck jplephem asdf

RUN conda install astropy
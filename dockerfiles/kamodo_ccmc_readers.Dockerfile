FROM ensemble/kamodo

# Install latest kamodo
ADD . /kamodo_ccmc

WORKDIR /kamodo_ccmc

RUN pip install -e .
# docker build -t asherp/kamodo -f API.Dockerfile .

FROM continuumio/miniconda3:latest
LABEL maintainer "Asher Pembroke <apembroke@predsci.com>"

RUN conda install python=3.7

# RUN conda install jupyter
RUN pip install antlr4-python3-runtime


# # Add Tini. Tini operates as a process subreaper for jupyter. This prevents kernel crashes.
# ENV TINI_VERSION v0.6.0
# ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
# RUN chmod +x /usr/bin/tini
# ENTRYPOINT ["/usr/bin/tini", "--"]

# need to pin this version for api
RUN pip install sympy==1.5.1

# Keep plotly at lower api
RUN pip install plotly==4.7.1

# Install latest kamodo
ADD . /kamodo

# RUN git clone https://github.com/asherp/kamodo.git
RUN pip install -e kamodo

RUN conda install jupyter

WORKDIR kamodo

# CMD ["kamodo-serve"]

CMD ["jupyter", "notebook", "./docs/notebooks", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]

#####
# For Jupyter notebook interaction, use:
#	docker run -p 8888:8888 dezeeuw/kamodo
# For command line interaction, use:
#	docker run -it dezeeuw/kamodo /bin/bash
#   -above, with current working directory mounted in container, use
#	docker run -it --mount type=bind,source="$(pwd)",destination=/local,consistency=cached  dezeeuw/kamodo /bin/bash
#   -above, with persistent disk space, use
#	docker run -it --mount source=kamododisk,target=/kdisk dezeeuw/kamodo /bin/bash
#
# Persistent disk space command
#	docker volume create kamododisk
#
# docker build . -t dezeeuw/kamodo 

FROM continuumio/miniconda3:latest
LABEL maintainer "Darren De Zeeuw <darrens@umich.edu>"

RUN conda install jupyter
RUN pip install antlr4-python3-runtime


# Add Tini. Tini operates as a process subreaper for jupyter. This prevents kernel crashes.
ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]

# Install latest kamodo
COPY . kamodo
RUN pip install -e kamodo

WORKDIR kamodo

CMD ["python", "kamodo/cli/api.py"]

# CMD ["jupyter", "notebook", "./docs/notebooks", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]

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

# Start from the fredericklab base container
FROM fredericklab/basecontainer:v0.4.8

# Installing precomputed python packages
RUN mamba install -y pillow 

# Install micscriptlib
COPY . /src/micscriptlib
RUN cd /src/micscriptlib && \
    pip install . && \
    rm -rf /src/micscriptlib/build /src/micscriptlib/dist

# clean up
RUN mamba clean -y --all
RUN pip cache purge

# Create a shared $HOME directory
RUN useradd -m -s /bin/bash -G users micscriptlib
WORKDIR /home/micscriptlib
ENV HOME="/home/micscriptlib"

ENV IS_DOCKER_8395080871=1

RUN ldconfig
WORKDIR /tmp/

# set a non-root user
USER micscriptlib

ARG VERSION
ARG BUILD_DATE
ARG VCS_REF

LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="micscriptlib" \
      org.label-schema.description="micscriptlib - a set of routines to simplify writing command line neuroimaging workflows" \
      org.label-schema.url="http://nirs-fmri.net" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/bbfrederick/micscriptlib" \
      org.label-schema.version=$VERSION 

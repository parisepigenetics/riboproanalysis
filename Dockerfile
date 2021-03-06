############################################################
# Dockerfile to build RiboProAnalysis container images
# Based on jpetazzo/dind
############################################################

#TODO this way is obsolete now... there is an official DockerInDocker image, we can pull it from DockerHub and to our job with that.

# Set the base image to Dind
FROM jpetazzo/dind

ENV http_proxy http://gwebproxy.univ-paris-diderot.fr:3128/
ENV https_proxy https://gwebproxy.univ-paris-diderot.fr:3128/

# File Author / Maintainer
MAINTAINER Alexandra Bomane

# Update repository
RUN apt-get update

# Install GNU parallel and matplotlib (Python library)
RUN apt-get install --yes parallel python-matplotlib

# Add application scripts
RUN mkdir /usr/local/riboproanalysis
ADD . /usr/local/riboproanalysis

# Mount the volume containing scripts for all created images
VOLUME /usr/local/riboproanalysis

# Create volumes to mount user's input files
RUN mkdir /rRNAindexdirectory
RUN mkdir /genomeindexdirectory

VOLUME /rRNAindexdirectory
VOLUME /genomeindexdirectory

RUN mkdir -p /usr/local/bin/riboproanalysis/
ADD . /usr/local/bin/riboproanalysis/

VOLUME /usr/local/bin/riboproanalysis/

# Update $PATH
ENV PATH $PATH:/usr/local/bin/riboproanalysis/RScripts/
ENV PATH $PATH:/usr/local/bin/riboproanalysis/PythonScripts/
ENV PATH $PATH:/usr/local/bin/riboproanalysis/

ENV http_proxy ""
ENV https_proxy ""

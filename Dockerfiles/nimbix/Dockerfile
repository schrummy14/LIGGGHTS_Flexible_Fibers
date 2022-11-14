FROM ubuntu:focal-20220531
LABEL maintainer="Matthewwschramm@gmail.com"

# Update SERIAL_NUMBER to force rebuild of all layers (don't use cached layers)
ARG SERIAL_NUMBER
ENV SERIAL_NUMBER ${SERIAL_NUMBER:-20200205.1000}

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata && apt-get install -y build-essential openmpi-bin libopenmpi-dev python-dev git bc paraview libvtk7-dev

# For general use, you may comment out the next RUN command
RUN apt-get -y install curl && curl -H 'Cache-Control: no-cache' \
        https://raw.githubusercontent.com/nimbix/image-common/ubuntu-focal-support/install-nimbix.sh \
        | bash -s -- --setup-nimbix-desktop --image-common-branch ubuntu-focal-support

# Build LIGGGHTS
RUN git clone https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers.git /opt/LIGGGHTS
RUN cd /opt/LIGGGHTS/src && make clean-auto && make -j$(python3 -c 'import multiprocessing as mp; print(int(mp.cpu_count() * 1.5))') auto
RUN ln -s /opt/LIGGGHTS/src/lmp_auto /usr/local/bin/liggghts_flexible_fibers

# Add PushToCompute Work Flow Metadata
ADD ./NAE/AppDef.json /etc/NAE/AppDef.json

# For general use, you do not need to expose these ports
# Expose port 22 for local JARVICE emulation in docker
EXPOSE 22

# for standalone use
EXPOSE 5901
EXPOSE 443

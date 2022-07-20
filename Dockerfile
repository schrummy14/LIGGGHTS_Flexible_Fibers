FROM ubuntu:focal-20220531
LABEL maintainer="Matt Schramm"

# Update SERIAL_NUMBER to force rebuild of all layers (don't use cached layers)
ARG SERIAL_NUMBER
ENV SERIAL_NUMBER ${SERIAL_NUMBER:-20200205.1000}

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata && apt-get install -y build-essential openmpi-bin libopenmpi-dev python-dev git bc paraview libvtk7-dev

RUN useradd --shell /bin/bash liggghts
RUN mkdir -p /home/liggghts/LIGGGHTS && chown -R liggghts:liggghts /home/liggghts
USER liggghts

RUN git clone https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers.git /home/liggghts/LIGGGHTS
RUN cd /home/liggghts/LIGGGHTS/src && pwd && make clean-auto && make -j$(python3 -c 'import multiprocessing as mp; print(int(mp.cpu_count() * 1.5))') auto
RUN mkdir -p /home/liggghts/.local/bin/
RUN ln -s /home/liggghts/LIGGGHTS/src/lmp_auto /home/liggghts/.local/bin/liggghts
RUN echo "export PATH=$PATH:$HOME/.local/bin" >> $HOME/.bashrc

WORKDIR /home/liggghts/

# Expose port 22 for local JARVICE emulation in docker
EXPOSE 22

# for standalone use
EXPOSE 5901
EXPOSE 443

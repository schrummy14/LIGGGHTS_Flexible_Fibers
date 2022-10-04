FROM ubuntu:22.04
LABEL maintainer="matthewwschramm@gmail.com"

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata && apt-get install -y build-essential cmake nano ninja-build openmpi-bin libopenmpi-dev python3-dev git bc paraview libvtk7-dev

# Add liggghts user
RUN useradd --shell /bin/bash liggghts
RUN mkdir -p /home/liggghts/LIGGGHTS && chown -R liggghts:liggghts /home/liggghts
USER liggghts

# Download LIGGGHTS
RUN git clone https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers.git /home/liggghts/LIGGGHTS

# Edit the utils file
RUN cd /home/liggghts/LIGGGHTS/src && sed 's/return static_cast< std::ostringstream & >(( std::ostringstream() << std::dec << a ) ).str();/return std::to_string(a);/' utils.h > notAfile.h && mv notAfile.h utils.h

# Build with cmake using ninja
RUN cd /home/liggghts/LIGGGHTS && mkdir BUILD && cd BUILD && cmake -GNinja ../src/. && ninja

# add local bin directory, add symbolic link, add directory to PATH
RUN mkdir -p /home/liggghts/.local/bin/
RUN ln -s ~/LIGGGHTS/BUILD/liggghts ~/.local/bin/.
RUN echo "export PATH=$PATH:$HOME/.local/bin" >> $HOME/.bashrc
RUN echo "export OMPI_MCA_btl_vader_single_copy_mechanism=none" >> $HOME/.bashrc

WORKDIR /home/liggghts/

# mpirun -n 8 liggghts -in in.chute_wear |& grep -v "Read -1"
# docker run -it -v $(pwd):/Data liggghts_22_04

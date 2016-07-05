#Example of docker configuration file
FROM ubuntu:14.04
RUN apt-get update && apt-get install -y \
     vim \
     htop \
     htop \
     libpng-dev \
     libjpeg-dev \
     libtiff-dev \
     libxxf86vm1 \
     libxxf86vm-dev \
     libxi-dev \
     libxrandr-dev \
     git \
     graphviz \
     cmake \
     g++ && \
     git clone --recursive https://github.com/openMVG/openMVG.git /tmp/openMVG && \
     cd /tmp/openMVG && \
     mkdir build && \
     cd build && \
     cmake ../src -DCMAKE_BUILD_TYPE=RELEASE -DOpenMVG_BUILD_EXAMPLES=ON -DOpenMVG_MAKE_TEST=ON \
     && make all \
     && make install \
 CMD ["make test"]

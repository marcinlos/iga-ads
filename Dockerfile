FROM ubuntu:16.04

ARG APP_DIR=/app
# Follows http://home.agh.edu.pl/~paszynsk/FastSmooth/short_course2.pdf

# Install common libraries
RUN apt update \
  && apt -y install wget build-essential cmake gfortran liblapack-dev libblas-dev libboost-all-dev gnuplot

# Install Galois
RUN wget http://iss.ices.utexas.edu/projects/galois/downloads/Galois-2.2.1.tar.gz \
   && tar xzvf Galois-*.tar.gz \
   && cd Galois-2.2.1/build \
   && mkdir release \
   && cd release \
   && cmake -DSKIP_COMPILE_APPS=ON ../.. \
   && make && make install

# Install libunittest
RUN wget https://datapacket.dl.sourceforge.net/project/libunittest/libunittest-9.3.5.tar.gz \
  && tar -xvf libunittest-*.tar.gz \
  && cd libunittest-*/ \
  && ./configure && make && make install

RUN mkdir -p $APP_DIR
COPY . $APP_DIR/
WORKDIR $APP_DIR
RUN cmake . && make
ENTRYPOINT ['cahn_hilliard']

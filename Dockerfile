FROM silkeh/clang:12 AS iga-ads-env

RUN apt-get update
RUN apt-get install -y libblas-dev
RUN apt-get install -y liblapack-dev
RUN apt-get install -y libboost-all-dev
RUN apt-get install -y git
RUN apt-get install -y curl

FROM iga-ads-env AS iga-ads-build

WORKDIR /src
COPY . .
RUN mkdir /build
RUN mkdir /ads
ENV DEPS=/ads/deps
# Make sure line endings are in unix style
RUN cat scripts/install-dependencies.sh | tr -d '\r' > scripts/install-dependencies.sh
RUN cat scripts/docker-start.sh | tr -d '\r' > /ads/docker-start.sh
RUN chmod +x /ads/docker-start.sh

RUN scripts/install-dependencies.sh /build/deps-build ${DEPS}
RUN cmake -S . -B /build \
    -D CMAKE_BUILD_TYPE=Release \
    -D ADS_USE_GALOIS=ON \
    -D ADS_USE_MUMPS=OFF \
    -D ADS_BUILD_TOOLS=OFF \
    -D CMAKE_PREFIX_PATH="${DEPS}"
RUN cmake --build /build --parallel --target ADS

RUN cmake --install /build

RUN rm -fr /build
RUN rm -fr /src

FROM iga-ads-build AS iga-ads

RUN echo "alias compile='cmake -S /src -B /build -D CMAKE_PREFIX_PATH=\"\${DEPS};/src/libtorch\" && make -C /build'" >> ~/.bashrc
WORKDIR /src
VOLUME /src
VOLUME /build

ENTRYPOINT [ "bash" , "--init-file", "/ads/docker-start.sh"]

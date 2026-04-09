FROM ubuntu:24.04

WORKDIR /code

ENV DEBIAN_FRONTEND=noninteractive

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        libzstd-dev \
        git \
        curl \
        ca-certificates \
        just \
        gfortran \
        g++ \
        ninja-build \
        liblapack-dev \
        libboost-all-dev \
        llvm-dev \
        libmumps-dev

# Install cmake
RUN --mount=type=bind,source=scripts/,target=scripts/ \
    scripts/install-cmake.sh

ENV CMAKE_GENERATOR=Ninja \
    CMAKE_COLOR_DIAGNOSTICS=ON

RUN --mount=type=bind,source=scripts/install-dependencies.sh,target=scripts/install-dependencies.sh,z \
    scripts/install-dependencies.sh /deps-build /deps

COPY . .

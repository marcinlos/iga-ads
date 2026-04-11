FROM ubuntu:25.10

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

# Install recent clang tools
RUN cat >> /etc/apt/sources.list <<-EOF
    deb http://apt.llvm.org/questing/ llvm-toolchain-questing-22 main
    deb-src http://apt.llvm.org/questing/ llvm-toolchain-questing-22 main
EOF

RUN curl -fLo /etc/apt/trusted.gpg.d/apt.llvm.org.asc https://apt.llvm.org/llvm-snapshot.gpg.key

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        clang-format-22 \
        clang-tidy-22 \
        && update-alternatives --install \
            /usr/bin/clang-format clang-format /usr/bin/clang-format-22 1 \
        && update-alternatives --install \
            /usr/bin/clang-tidy clang-tidy /usr/bin/clang-tidy-22 1

# Install prek
COPY --from=ghcr.io/j178/prek:v0.3.8 /prek /usr/local/bin/

# Install cmake
RUN --mount=type=bind,source=scripts/,target=scripts/ \
    scripts/install-cmake.sh

ENV CMAKE_GENERATOR=Ninja \
    CMAKE_COLOR_DIAGNOSTICS=ON

RUN --mount=type=bind,source=scripts/install-dependencies.sh,target=scripts/install-dependencies.sh \
    scripts/install-dependencies.sh /deps-build /deps

COPY . .

# Workaround for docker
RUN git config --global --add safe.directory /code

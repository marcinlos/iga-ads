FROM ubuntu:25.10

WORKDIR /code

ENV DEBIAN_FRONTEND=noninteractive

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        pkg-config \
        zip \
        unzip \
        git \
        curl \
        ca-certificates \
        just \
        gfortran \
        g++ \
        ninja-build \
        liblapack-dev \
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

# Setup vcpkg
RUN git clone https://github.com/microsoft/vcpkg.git --depth=1 /opt/vcpkg && \
    /opt/vcpkg/bootstrap-vcpkg.sh

ENV VCPKG_ROOT=/opt/vcpkg \
    VCPKG_FORCE_SYSTEM_BINARIES=ON

RUN ln -s "${VCPKG_ROOT}/vcpkg" /usr/local/bin/vcpkg

ENV CMAKE_GENERATOR=Ninja \
    CMAKE_COLOR_DIAGNOSTICS=ON \
    CMAKE_INSTALL_PREFIX=/opt/ads

# Populate vcpkg binary cache
RUN --mount=type=bind,source=vcpkg.json,target=vcpkg.json \
    --mount=type=bind,source=vcpkg-configuration.json,target=vcpkg-configuration.json \
    --mount=type=bind,source=vcpkg,target=vcpkg \
    vcpkg install && rm -rf vcpkg_installed

# Workaround for docker
RUN git config --global --add safe.directory /code

# Pre-install prek hooks
RUN --mount=type=bind,source=prek.toml,target=prek.toml \
    # prek needs a git repo \
    git init -b main . && \
    prek prepare-hooks && \
    rm -rf .git

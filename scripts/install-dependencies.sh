#!/usr/bin/env bash

read -r -d '' usage <<-EOF
usage: $(basename "$0") build_dir install_dir

Download, build and install ADS dependencies (fmt, Galois, Catch2)

Positional arguments:
 build_dir       directory where the sources are downloaded and built
 install_dir     installation prefix
EOF

if [[ $# -ne 2 ]]; then
    echo "${usage}"
    exit 1
fi

set -ex

IGA_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." &> /dev/null && pwd)

BUILD_DIR=${1}
INSTALL_DIR=${2}

BUILD_TYPE=Release

if [[ -d "${INSTALL_DIR}" ]]; then
  echo "Installation dir ${INSTALL_DIR} exists, skipping installing dependencies"
  exit
fi

mkdir -p "${BUILD_DIR}"
mkdir -p "${INSTALL_DIR}"

CMAKE_BUILD_PARALLEL_LEVEL=$(nproc)
export CMAKE_BUILD_PARALLEL_LEVEL


cd "${BUILD_DIR}"

# Install Galois
GALOIS_VER=6.0
git clone --branch release-${GALOIS_VER} --depth=1 --quiet https://github.com/IntelligentSoftwareSystems/Galois

git -C Galois apply "${IGA_DIR}/scripts/galois.patch"

mkdir -p Galois/build

cmake \
  -S Galois \
  -B Galois/build \
  -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -D CMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -D CMAKE_PREFIX_PATH="${INSTALL_DIR}" \
  -D GALOIS_ENABLE_DIST=OFF \
  -D BUILD_TESTING=OFF

cmake --build Galois/build --target galois_shmem
cmake --install Galois/build --component dev
cmake --install Galois/build --component lib

# Needed because Galois imported targets reference these
cmake --build Galois/build --target graph-convert
cmake --build Galois/build --target graph-convert-huge
cmake --install Galois/build --component tools

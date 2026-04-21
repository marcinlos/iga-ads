#! /usr/bin/bash

REPO="https://github.com/Kitware/CMake"
VER=4.3.2
CMAKE_URL="${REPO}/releases/download/v${VER}/cmake-${VER}-linux-x86_64.sh"
TARGET=$(mktemp /tmp/cmake-installer.XXXXX)

curl -fLo "${TARGET}" "${CMAKE_URL}"
bash "${TARGET}" --prefix=/usr/local --skip-license
rm "${TARGET}"

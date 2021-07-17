#!/usr/bin/env bash

set -eux

IFS=- read -r COMPILER VERSION <<< "${1,,}"

echo "Compiler: ${COMPILER}, version ${VERSION}"

case ${COMPILER} in
  gcc)
    cxx="g++-${VERSION}"
    # Matching version of Fortran standard library is needed
    sudo apt-get install "${cxx}" "libgfortran-${VERSION}-dev"
    echo "CXX=${cxx}" >> "${GITHUB_ENV}"
    ;;

  clang)
    # Prior to 7, packages have -major.minor suffix
    if (( VERSION < 7 )); then
      VERSION="${VERSION}.0"
    fi

    sudo apt-get install "clang-${VERSION}"
    echo "CXX=clang++-${VERSION}" >> "${GITHUB_ENV}"
    ;;

  *)
    echo "Error: unknown compiler ${COMPILER}"
    exit 1
    ;;
esac

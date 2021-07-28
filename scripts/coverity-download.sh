#!/usr/bin/env bash

set -eux

# Download Coverity Scan Tool
wget -q https://scan.coverity.com/download/cxx/linux64 \
    --post-data "token=${TOKEN}&project=marcinlos/iga-ads" \
    -O coverity.tgz

# Extract into a directory with simple name
DIR=coverity

mkdir -p "${DIR}"

tar xzf coverity.tgz \
    --strip-components 1 \
    --directory "${DIR}"

# Ensure cov-build is on $PATH
echo "PATH=$(readlink -f ${DIR}/bin):${PATH}" >> "${GITHUB_ENV}"

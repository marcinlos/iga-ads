#!/usr/bin/env bash

FILE="scan-data.tgz"

tar czf "${FILE}" build/cov-init

curl \
    --form project=marcinlos/iga-ads \
    --form token=${TOKEN} \
    --form email=marcin.los.91@gmail.com \
    --form file=${FILE} \
    --form version=develop \
    --form description="develop build" \
    https://scan.coverity.com/builds

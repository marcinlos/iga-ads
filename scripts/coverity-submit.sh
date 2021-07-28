#!/usr/bin/env bash

set -eux

FILE="scan-data.tgz"

tar czf "${FILE}" cov-int

curl \
    --form token=${TOKEN} \
    --form email=marcin.los.91@gmail.com \
    --form file=@${FILE} \
    --form version=${GITHUB_SHA} \
    --form description="develop build" \
    https://scan.coverity.com/builds?project=marcinlos%2Figa-ads

#!/usr/bin/env bash

# Script for parsing output of run-clang-tidy and counting the number of
# occurrences of individual warnings.

file=${1:-warnings}
root=$( git rev-parse --show-toplevel )

script="
    \~^${root}~!d
    s/.\+\[\([a-zA-Z.-]\+\)\]/\1/
    \~^${root}~d
"

sed "${script}" "${file}" | sort | uniq --count | sort -n --reverse

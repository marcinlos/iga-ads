#!/bin/bash

alias compile='cmake -S /src -B /build -D CMAKE_PREFIX_PATH="${DEPS}" && make -C /build'
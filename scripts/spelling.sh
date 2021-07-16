#!/usr/bin/env bash

SKIP='.git,external,*build*'
ALLOWED='dof'
MASK=7

codespell -q ${MASK} -S ${SKIP} -L ${ALLOWED}

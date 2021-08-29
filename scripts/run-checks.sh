#!/usr/bin/env bash

DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

"${DIR}/spelling.sh"
"${DIR}/copyright.py" check
"${DIR}/guards.py" check
"${DIR}/format.py" check

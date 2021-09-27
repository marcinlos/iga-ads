#!/usr/bin/env bash

# Arguments: <accuracy level> <bx> <by> [ADS build dir]
# Without the last argument, value of ADS_BUILD env variable is used

level="$1"
bx="$2"
by="$3"

n_base=16
p=2
c=1
P=2
C=0

n_coarse=$(python3 -c "print(int(${n_base} * 2 ** (${level} - 1)))")
n_fine=$(python3 -c "print(int(${n_base} * 2 ** ${level}))")

export ADS_BUILD="${ADS_BUILD:-$4}"

dir=$(mktemp --directory)

cd "${dir}" || exit

"${ADS_BUILD}/examples/erikkson_inverse" "${n_coarse}" "$p" "$c" "$P" "$C" "$bx" "$by" &> /dev/null
mv coeffs.data coarse.data

"${ADS_BUILD}/examples/erikkson_inverse" "${n_fine}" "$p" "$c" "$P" "$C" "$bx" "$by" &> /dev/null
mv coeffs.data fine.data

"${ADS_BUILD}/examples/inverse_postprocess" \
    coarse.data "${n_coarse}" "$p" "$c" \
    fine.data "${n_fine}" "$p" "$c"

rm -rf "${dir}"

#!/usr/bin/env bash

mkdir -p data-raw/logs
PROG="${1}"
RESULTS="${2:-nonlinear}"

run_test() {
    N="${1}"
    ELEMS="${2}"
    SCHEME="${3}"

    DT=$(python3 -c "print(1/${N})")
    OUTPUT="data-raw/logs/results-${ELEMS}-${SCHEME}-${N}.txt"

    rm -f -- *.data
    stdbuf --output=L \
        "${PROG}" \
            --elems "${ELEMS}" -p 2 \
            --scheme "${SCHEME}" \
            --threads 24 \
            --steps "${N}" \
            --dt "${DT}" \
        > "${OUTPUT}"

    tail -1 "${OUTPUT}"
    tail -1 "${OUTPUT}" >> "data-raw/${RESULTS}-${ELEMS}-${SCHEME}"
}

# for N in 1 2 4 8 16 32 64 128 256 512 1024 2048 4096; do
#     for ELEMS in 50 100 200; do
#         for SCHEME in PR strang-CN FE; do
#             echo "${N} steps ${SCHEME} n=${ELEMS}"
#             run_test "${N}" "${ELEMS}" "${SCHEME}";
#         done
#     done
# done
for N in 1 2 4 8 16 32 64 128 256 512 1024; do
    for ELEMS in 100; do
        for SCHEME in PR; do
            # echo "${N} steps ${SCHEME} n=${ELEMS}"
            run_test "${N}" "${ELEMS}" "${SCHEME}";
        done
    done
done

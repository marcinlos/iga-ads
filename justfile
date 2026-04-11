# Task definitions for the just command runner
# https://just.systems/man/en/

image := "iga-ads:latest"
tool := env("CONTAINER_TOOL", "docker")
build_dir := "/build"

@_default:
    just --list
    echo
    echo "Using {{tool}}"

# Build the development container
@image:
    {{tool}} build \
        -f Containerfile \
        -t {{image}} \
        .

# Start the development container
@shell:
    {{tool}} run \
        --rm \
        --interactive \
        --tty \
        --volume .:/code:z \
        {{image}} \
        bash

@config:
    cmake \
        --fresh \
        -S /code \
        -B {{build_dir}} \
        -D CMAKE_BUILD_TYPE=Release \
        -D ADS_USE_GALOIS=ON \
        -D ADS_USE_MUMPS=ON \
        -D CMAKE_PREFIX_PATH=/deps

@build CORES="$(nproc)":
    cmake \
        --build {{build_dir}} \
        --parallel {{CORES}}

# Run all the linter tools
@lint:
    prek run --all-files

# Check that public headers are self-sufficient
@verify-headers:
    cmake \
        --build {{build_dir}} \
        --target all_verify_interface_header_sets

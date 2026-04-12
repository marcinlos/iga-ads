# Task definitions for the just command runner
# https://just.systems/man/en/

image := "iga-ads:latest"
container_name := "iga-ads-dev"
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
@start: && shell
    {{tool}} run \
        --rm \
        --detach \
        --tty \
        --volume .:/code:z \
        --name {{container_name}} \
        {{image}}

# Stop the development container
@stop:
    {{tool}} kill {{container_name}}

# Start a shell in the development container
@shell:
    {{tool}} exec \
        --interactive \
        --tty \
        {{container_name}} \
        bash

@config:
    cmake \
        --fresh \
        -S /code \
        -B {{build_dir}} \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_TOOLCHAIN_FILE=etc/gcc-toolchain.cmake \
        -D ADS_USE_GALOIS=ON \
        -D ADS_USE_MUMPS=ON \
        -D CMAKE_PREFIX_PATH=/deps \
        -D CMAKE_INSTALL_LIBDIR=lib \
        -D CMAKE_INSTALL_PREFIX=/opt/ads

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

# Install git hooks running inside the development container
@install-hooks:
    {{tool}} exec \
        {{container_name}} \
        prek install
    mv .git/hooks/pre-commit .git/hooks/invoke-prek
    ln --symbolic --force \
        ../../hooks/pre-commit \
        .git/hooks/

# Remove git hooks installed by `install-hooks`
@clear-hooks:
    rm -f .git/hooks/{pre-commit,pre-commit.legacy,invoke-prek}

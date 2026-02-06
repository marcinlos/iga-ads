@image:
    podman build -t iga-ads:latest .

@shell:
    podman run \
        --rm \
        --interactive \
        --tty \
        --volume .:/code:z \
        iga-ads:latest \
        bash

@config:
    cmake \
        -S /code \
        -B /build \
        -D CMAKE_BUILD_TYPE=Release \
        -D ADS_USE_GALOIS=ON \
        -D ADS_USE_MUMPS=ON \
        -D CMAKE_PREFIX_PATH=/deps

@build:
    cmake --build /build -j $(nproc)

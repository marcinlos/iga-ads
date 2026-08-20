# Run iga-ads in docker

Requirements:
* [Docker](https://www.docker.com)

## 1. Build docker image
From git root directory run:

```
docker build --tag iga-ads .
```

## 2. Create your project
This directory contains simple project that can be built and run using docker image created in the previous step. 
When creating new project you can copy `CMakeLists.txt` file and adjust its content.
Moreover, prepare two directories:
* src:   with simualtion code and `CMakeLists.txt`
* build: for build artifacts (may be empty)

## 3. Run docker container

To build and run this example create container using command (from this directory):

```
docker run -it --rm -v $PWD/:/src/ -v $PWD/build/:/build/ iga-ads
```

To build and run your project with prepared directories `/project/src` and `/project/build` run:

```
docker run -it --rm -v /project/src/:/src/ -v /project/build/:/build/ iga-ads
```

## 4. Compile
From docker bash run:

```
compile
```

## 5. Run
Compiled program is written in `/build` with name specified in `CMakeLists.txt`. To run example from this directory:

```
cd /build
./heat_1d
```
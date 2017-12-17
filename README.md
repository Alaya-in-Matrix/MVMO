# README

## About

This is a c++ implementation of the Mean Variance Mapping Optimization
algorithm(MVMO). The algorithm has pretty good performance in CEC14 expensive
black-box single-objective optimization competition.

## Dependency

- CMake for compiling 
- Eigen for matrix operations
- Nlopt for local search

## Install 

```bash
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/of/installation -DEigen3_DIR=/path/of/Eigen/share/eigen3/cmake -DNLOPT_PATH=/path/to/nlopt
make install
```

## DEMO

- TODO

## Origin algorithm

- The [MVMO home page](https://www.uni-due.de/mvmo/)
- [Mean Variance Mapping Optimization](https://www.uni-due.de/imperia/md/content/mvmo/background_mvmo.pdf)

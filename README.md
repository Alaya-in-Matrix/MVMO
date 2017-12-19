# README

## About

This is a c++ implementation of the Mean Variance Mapping Optimization
algorithm(MVMO). The algorithm has pretty good performance in CEC14 expensive
black-box single-objective optimization competition.

## Dependency

- CMake for compiling 
- Eigen for matrix operations

## Install 

```bash
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/of/installation -DEigen3_DIR=/path/of/Eigen/share/eigen3/cmake
make install
```

## DEMO

- See demo.cpp

## Origin algorithm

- The [MVMO home page](https://www.uni-due.de/mvmo/)
- [Mean Variance Mapping Optimization](https://www.uni-due.de/imperia/md/content/mvmo/background_mvmo.pdf)
- Erlich, István, et al. "Solving the IEEE-CEC 2014 expensive optimization test
  problems by using single-particle MVMO." Evolutionary Computation (CEC), 2014
  IEEE Congress on. IEEE, 2014.

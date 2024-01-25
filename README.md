# A Smooth Backfitting (SBF) header file library written in C++

This is a C++ implementation of various smooth backfitting techniques, see e.g. ["The Existence and Asymptotic properties of a Backfitting Projection Algorithm under Weak Conditions"](https://projecteuclid.org/journals/annals-of-statistics/volume-27/issue-5/The-existence-and-asymptotic-properties-of-a-backfitting-projection-algorithm/10.1214/aos/1017939138.pdf)

## Installation

Simply add the header files from the `/include` folder into your own directory.

The only dependency is the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library for matrix operations.

## When you should (and shouldn't) use this implementation

To the best of my knowledge, the only other public implementation of the smooth backfitting algorithm is in the R package [`fdapace`](https://www.rdocumentation.org/packages/pcalg/versions/2.7-1). 

Thus, **this implementation might be for you if**:

- you want an easy-to-use and dependency-light implementation of the smooth backfitting procedure, or
- you want a computationally fast implementation of the smooth backfitting method or variations thereof
  
**You should not use this implementation if:**

- you want to call these methods from Python. In this case, please use the wrapper library [`SmoothBackfittingPy`](https://github.com/AsgerMorville/SmoothBackfittingPy) which calls these C++ methods using an extension module

## Running the algorithm

Here is a simple example to get you started. Observations from an additive function with a 2-dimensional domain are fitted and the fitted values are printed.

```c++
#include <iostream>
#include <cmath>
#include <random>

#include "Eigen/Dense"
#include "smooth_backfitting_core.h"
#include "additive_function.h"
#include "utils/random_generator.h"

typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXXd Array;
typedef Eigen::MatrixXd Matrix;


int main(){
    // Set random number generation seed
    std::mt19937 mt(1);

    // Set size parameters
    size_t n = 20; // number of data points
    size_t d = 2;  // dimension

    // Generate design matrix X
    Matrix X = uniformMatrix(n, d, mt);

    // Generate observations
    Vector noise = stdNormalVector(n, mt);
    Vector Y = Vector::Zero(n);
    for (int i = 0; i < n; i++){
        Y(i) = std::cos(X(i,0)) + std::sin(X(i,1)) + noise(i);
    }

    // Obtain smooth backfitting estimates
    AddFunction addFunc = SBF(Y,X);

    // Use the .predict method to evaluate the estimated additive function. For example,
    std::cout << "Fitted points: \n" << addFunc.predict(X) << "\n";

    return 0;
};
```


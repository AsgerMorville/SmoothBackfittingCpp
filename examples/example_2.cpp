// This file is part of SmoothBackfittingCpp, a header file library for smooth backfitting methods in C++
//
// Copyright (C) 2023-2024 <asgermorville@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <cmath>
#include <random>

#include "Eigen/Dense"
#include "partially_linear_SBF.h"
#include "additive_function.h"
#include "utils/random_generator.h"

typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXXd Array;
typedef Eigen::MatrixXd Matrix;

double f_1(double x){
    return std::sin(x);
}

double f_2(double x){
    return std::cos(x);
}


int main(){
    // Example of the partially linear SBF method

    // Set random number generation seed
    std::mt19937 mt(1);

    // Set size parameters
    size_t n = 20; // number of data points
    size_t d_x = 2; // dimension of the linear part
    size_t d_z = 2;  // dimension of the non-parametric part

    // Generate design matrices X and Z
    Matrix X = bernoulliMatrix(n, d_x, mt);
    Matrix Z = uniformMatrix(n, d_z, mt);

    // Generate observations
    Vector noise = stdNormalVector(n, mt);
    Vector Y = Vector::Zero(n);
    for (int i = 0; i < n; i++){
        Y(i) = X(i,0)*3 + X(i,1)*(-2) + f_1(Z(i,0)) + f_2(Z(i,1)) + noise(i);
    }

    // Obtain smooth backfitting estimates as a partially linear additive function object
    PartAddFunction pAddFunc = PL_SBF(Y, X, Z);

    // Use the .predict method to evaluate estimated additive function. For example,
    std::cout << "Fitted points: \n" << pAddFunc.predict(X, Z) << "\n";

    return 0;
};
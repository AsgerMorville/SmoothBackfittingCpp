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
#include "smooth_backfitting_core.h"
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
    // Additive function example

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
        Y(i) = f_1(X(i,0)) + f_2(X(i,1)) + noise(i);
    }

    // Obtain smooth backfitting estimates
    AddFunction addFunc = SBF(Y,X);

    // Use the .predict method to evaluate estimated additive function. For example,
    std::cout << "Fitted points: \n" << addFunc.predict(X) << "\n";

    return 0;
};

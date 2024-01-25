// This file is part of SmoothBackfittingCpp, a header file library for smooth backfitting methods in C++
//
// Copyright (C) 2023-2024 <asgermorville@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SMOOTH_BACKFITTING_LIBRARY_RANDOM_GENERATOR_H
#define SMOOTH_BACKFITTING_LIBRARY_RANDOM_GENERATOR_H

#include <Eigen/Dense>
#include <random>

typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXXd Array;
typedef Eigen::MatrixXd Matrix;

Matrix uniformMatrix(size_t n, size_t d, std::mt19937 seed){
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    Matrix X = Matrix::Zero(n,2);
    for (int j = 0; j < d; j++){
        for (int i = 0; i < n; i++){
            X(i,j) = unif(seed);
        }
    }
    return X;
}

Matrix bernoulliMatrix(size_t n, size_t d, std::mt19937 seed){
    std::bernoulli_distribution bern(0.5);
    Matrix X = Matrix::Zero(n,2);
    for (int j = 0; j < d; j++){
        for (int i = 0; i < n; i++){
            X(i,j) = bern(seed);
        }
    }
    return X;
}


Vector stdNormalVector(size_t n, std::mt19937 seed){
    std::normal_distribution<double> norm(0,1);
    Vector Y = Vector::Zero(n);
    for (int i = 0; i < n; i++){
        Y(i) = norm(seed);
    }
    return Y;
}


#endif //SMOOTH_BACKFITTING_LIBRARY_RANDOM_GENERATOR_H

//
// Created by Asger Morville on 2024/01/07.
//
#include <Eigen/Dense>
#include <random>
#include <iostream>

int main() {
    // Dimensions of the matrix
    int rows = 10;
    int cols = 5;

    // Random number generation setup
    std::random_device rd;  // Seed for the random number engine
    std::default_random_engine eng(rd()); // Random number engine
    std::uniform_int_distribution<> distr(0, 1); // Distribution for 0 or 1

    // Create and initialize the matrix with random 0 or 1
    Eigen::MatrixXd mat = Eigen::MatrixXd::NullaryExpr(rows, cols, [&]() { return distr(eng); });

    // Print the matrix
    std::cout << mat << std::endl;

    return 0;
}
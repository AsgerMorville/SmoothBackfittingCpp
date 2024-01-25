//
// Created by Asger Morville on 2024/01/12.
//

#include "partially_linear_SBF.h"
#include <random>
typedef Eigen::ArrayXXd Array;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

int main() {
    int n = 100;
    int d = 5;

    // Random number generation setup
    std::random_device rd;  // Seed for the random number engine
    std::default_random_engine eng(rd()); // Random number engine
    std::uniform_int_distribution<> distr(0, 1); // Distribution for 0 or 1

    // Create and initialize the matrix with random 0 or 1
    Matrix X = Matrix::NullaryExpr(n, d, [&]() { return distr(eng); });
    Matrix Z = Matrix::Random(n, d);
    Vector Y = Vector::Random(n);


    /*
    auto t0 = std::chrono::steady_clock::now();
    Vector fitted = PL_SBF(Y, X, Z);
    auto t1 = std::chrono::steady_clock::now();
    //std::cout << add_func_test.predict(X) << "\n";
    std::chrono::duration<double, std::milli> timee = t1 - t0;
    std::cout << (timee / 1000.0).count() << " seconds" << "\n";
    return 0;
    */
}
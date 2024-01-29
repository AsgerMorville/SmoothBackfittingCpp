//
// Created by Asger Morville on 2024/01/07.
//
#include <iostream>
#include "smooth_backfitting_core.h"
#include <Eigen/Dense>
#include <chrono>
typedef Eigen::ArrayXXd Array;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

int main(){
    /*
    int n = 2000;
    int d = 20;
    Matrix X = Matrix::Random(n,d);
    Vector Y = Vector::Random(n);

    auto t0 = std::chrono::steady_clock::now();
    AddFunction add_func_test = SBF(Y,X);
    auto t1 = std::chrono::steady_clock::now();

    //std::cout << add_func_test.predict(X) << "\n";
    std::chrono::duration<double, std::milli> timee = t1-t0;
    std::cout << (timee/1000.0).count() << " seconds"  << "\n";
    */
    int n = 10;
    int d = 2;
    Matrix X(n,d);
    Array m(n, d);
    Vector input_point(2);
    Array input_array(3,d);
    Vector Y(n);

    Y << 0.6439, -1.0304, -0.1212, -1.2307, -0.8734, -0.1489, -1.0219,
            0.1086, -0.8173, -1.3086;
    X.col(0) << -0.7663,  0.4802, -0.6479,  1.3356, -0.8366,  0.1365,  0.4252, -0.4268,  0.55  ,  0.281;
    X.col(1) << -0.1501, -0.7642,  0.0367,  0.5074, -0.9581,  0.9965,  0.7849, 0.4374,  0.7842, -0.5952;


    AddFunction add_func_test = SBF(Y,X);

    std::cout << "Y: " << Y << "\n";
    std::cout << "X : " << X << "\n";
    std::cout << "Xi evaluations: \n" << add_func_test.predict(X) << "\n";
    // This should print: [ 0.3038 -1.3868  0.2132 -1.2521 -0.7052 -0.4725 -0.7166  0.027  -0.764
    // -1.2051]


    /*
    AddFunction test = SBF(Y,X);
    Array input_point2(1,d);
    input_point2 << 0.2, 0.5;
    std::cout << "X_i evals: " << test.predict(X) << "\n";
    //std::cout << "x evals: " << test.xEvaluations() << "\n";
    //std::cout << X*X.transpose();
    //Array sdtest(4,2);
    //sdtest << 0.2856,  0.8851, -0.7544,  1.2529, 0.5129, -0.2981, 0.4885, -0.0756;
    //std::cout << "input \n" <<  sdtest << "\n";
    //std::cout << "hinit test \n" <<  hInitialize(sdtest) << "\n";

    //sbfFitter testfitt;
    //AddFunction test2 = testfitt.fit(Y,X);
    //std::cout << "objectification: " << test2.predict(X) << "\n";

    Vector output = Vector::Zero(n);
    Array X2 = X.transpose();
    std::cout << "Before FITTING: " << output << "\n";
    sbfWrapper(Y.data(), X2.data(),output.data(), n, d);
    std::cout << "AFTER FITTING: " << output << "\n";
    */
    return 0;
};
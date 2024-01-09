//
// Created by Asger Morville on 2024/01/07.
//
#include <iostream>
#include <Eigen/Dense>
#include "smooth_backfitting_core.h"
typedef Eigen::ArrayXXd Array;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

int main(){
    int n = 10;
    int n0 = 5;
    int d = 2;
    Array X(n,d);
    Array Xi(n0, d);
    Array m(n, d);
    Vector input_point(2);
    Array input_array(3,d);
    Vector Y(n);
    double y_mean = 2.3;
    X << 0.417, 0.7203, 0.0001, 0.3023, 0.1468, 0.0923, 0.1863, 0.3456, 0.3968, 0.5388, 0.4192, 0.6852,
            0.2045, 0.8781, 0.0274, 0.6705, 0.4173, 0.5587, 0.1404, 0.1981;
    Xi << 0.8007, 0.9683, 0.3134, 0.6923, 0.8764, 0.8946, 0.085 , 0.0391, 0.1698, 0.8781;
    m << -0.1724, -0.8779, 0.0422,  0.5828, -1.1006,  1.1447, 0.9016,  0.5025, 0.9009, -0.6837,
            -0.1229, -0.9358, -0.2679,  0.5304, -0.6917, -0.3968, -0.6872, -0.8452, -0.6712, -0.0127;
    input_point << 0.2, 0.5;
    input_array << 0.5741, 0.1467, 0.5893, 0.6998, 0.1023, 0.4141;

    Y << 2.2683, 1.7119, 1.1314, 1.8573, 2.0992, 2.0611, 2.6933, 2.2075, 2.3477, 1.3244;
    /*
    AddFunction add_func_test = SBF(Y,X);

    std::cout << "TESTER eval: " << testerr << "\n";
    std::cout << "TRUE eval: " << 1.3484 <<"\n";

    std::cout << "TESTER predict: " << tester_array << "\n";
    std::cout << "TRUE predict: " << "2.1644, 1.2413, 2.1312" << "\n";
     */

    AddFunction test = SBF(Y,X);
    Array input_point2(1,d);
    input_point2 << 0.2, 0.5;
    std::cout << "predict: " << test.predict(input_point2) << "\n";
    std::cout << "x evals: " << test.xEvaluations() << "\n";
    //std::cout << X*X.transpose();
    //Array sdtest(4,2);
    //sdtest << 0.2856,  0.8851, -0.7544,  1.2529, 0.5129, -0.2981, 0.4885, -0.0756;
    //std::cout << "input \n" <<  sdtest << "\n";
    //std::cout << "hinit test \n" <<  hInitialize(sdtest) << "\n";

    sbfFitter testfitt;
    AddFunction test2 = testfitt.fit(Y,X);
    std::cout << "objectification: " << test2.xEvaluations() << "\n";
    return 0;
};
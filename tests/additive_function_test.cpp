//
// Created by Asger Morville on 2024/01/07.
//
#include <iostream>
#include <Eigen/Dense>
#include "additive_function.h"
typedef Eigen::ArrayXXd Array;
typedef Eigen::VectorXd Vector;

int main(){
    int n = 10;
    int n0 = 5;
    int d = 2;
    Array X(n,d);
    Array Xi(n0, d);
    Array m(n, d);
    Vector input_point(2);
    Array input_array(3,d);
    double y_mean = 2.3;
    X << 0.417, 0.7203, 0.0001, 0.3023, 0.1468, 0.0923, 0.1863, 0.3456, 0.3968, 0.5388, 0.4192, 0.6852,
        0.2045, 0.8781, 0.0274, 0.6705, 0.4173, 0.5587, 0.1404, 0.1981;
    Xi << 0.8007, 0.9683, 0.3134, 0.6923, 0.8764, 0.8946, 0.085 , 0.0391, 0.1698, 0.8781;
    m << -0.1724, -0.8779, 0.0422,  0.5828, -1.1006,  1.1447, 0.9016,  0.5025, 0.9009, -0.6837,
        -0.1229, -0.9358, -0.2679,  0.5304, -0.6917, -0.3968, -0.6872, -0.8452, -0.6712, -0.0127;
    input_point << 0.2, 0.5;
    input_array << 0.5741, 0.1467, 0.5893, 0.6998, 0.1023, 0.4141;


    AddFunction test_func = AddFunction(X,m,Xi,y_mean);
    test_func.get_m_points();
    double testerr = test_func.eval(input_point);
    Array tester_array = test_func.predict(input_array);

    std::cout << "TESTER eval: " << testerr << "\n";
    std::cout << "TRUE eval: " << 1.3484 <<"\n";

    std::cout << "TESTER predict: " << tester_array << "\n";
    std::cout << "TRUE predict: " << "2.1644, 1.2413, 2.1312" << "\n";
    return 0;
};
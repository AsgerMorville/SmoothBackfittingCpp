//
// Created by Asger Morville on 2024/01/07.
//
#include <iostream>
#include <Eigen/Dense>
#include "smooth_backfitting_core.h"
#include "additive_function.h"

typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXXd Array;

int main(){
    // Additive function example

    size_t n = 20;
    size_t d = 1;

    Vector xx = Vector::LinSpaced(n,0,1);
    Array x_pp = xx.replicate(1,d);
    Array m_pp = xx.replicate(1,d);
    Array X_p = Vector::Random(n);
    double y_meann = 2.33;

    AddFunction tester = AddFunction(x_pp, m_pp, X_p, y_meann);
    //Array x_p, Array m_p, Array X_i, double y_mean
    std::cout << "X_p: \n" << X_p << "\n";
    std::cout << "tester: " << tester.get_m_points() << "\n";
    return 0;
};

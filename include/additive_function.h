//
// Created by Asger Morville on 2024/01/07.
//

#ifndef SMOOTH_BACKFITTING_LIBRARY_ADDITIVE_FUNCTION_H
#define SMOOTH_BACKFITTING_LIBRARY_ADDITIVE_FUNCTION_H

#include <Eigen/Dense>

typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXXd Array;

class AddFunction{
private:
    Array m_x_points;
    Array m_m_points;
    size_t m_d;
    size_t m_m;
    double m_y_mean;
    Array m_X_i;
public:
    AddFunction(Array x_p, Array m_p, Array X_i, double y_mean)
            : m_x_points(std::move(x_p)), m_m_points(std::move(m_p)), m_X_i(std::move(X_i)), m_y_mean(y_mean) {
        m_d = m_x_points.cols();
        m_m = m_x_points.rows();
    }
    Array get_m_points(){
        return m_m_points;
    };
    //double eval_function(Vector x_points);
    //Vector predict(Array input);
    //Vector X_i_evaluations;
};

#endif //SMOOTH_BACKFITTING_LIBRARY_ADDITIVE_FUNCTION_H

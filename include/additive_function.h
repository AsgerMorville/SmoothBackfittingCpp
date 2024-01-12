//
// Created by Asger Morville on 2024/01/07.
//

#ifndef SMOOTH_BACKFITTING_LIBRARY_ADDITIVE_FUNCTION_H
#define SMOOTH_BACKFITTING_LIBRARY_ADDITIVE_FUNCTION_H

#include <Eigen/Dense>

typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXXd Array;
typedef Eigen::MatrixXd Matrix;

class AddFunction{
private:
    Array m_x_points;
    Array m_m_points;
    size_t m_d;
    size_t m_m;
    double m_y_mean;
public:
    AddFunction(Array x_p, Array m_p, double y_mean)
            : m_x_points(std::move(x_p)), m_m_points(std::move(m_p)), m_y_mean(y_mean) {
        m_d = m_x_points.cols();
        m_m = m_x_points.rows();
    }
    Array get_m_points(){
        std::cout << "mdd: " << m_d << "\n";
        return m_m_points;
    };

    /**
     * Evaluates the d-dimensional additive function at the input point
     * @param x-coordinates of the location where f is to be evaluated
     * @return value of the additive function at x
     */
    double eval(Vector input){
        double output = 0;
        int picked_index;
        for (int j =0; j < m_d; j++){
            // Find minimizing index for each dimension
            Vector difference_vector = (m_x_points.col(j).array() - input(j)).abs();
            difference_vector.minCoeff(&picked_index);
            output += m_m_points(picked_index,j);
        }
        return output+m_y_mean;
    }

    Vector predict(Array &input){
        int input_length = input.rows();
        Vector output(input_length);
        for (int i = 0; i < input_length; i++){
            output[i] = eval(input.row(i));
        }
        return output;
    }
    Vector predict(Matrix &input2){
        Array input = input2.array();
        int input_length = input.rows();
        Vector output(input_length);
        for (int i = 0; i < input_length; i++){
            output[i] = eval(input.row(i));
        }
        return output;
    }
};

#endif //SMOOTH_BACKFITTING_LIBRARY_ADDITIVE_FUNCTION_H

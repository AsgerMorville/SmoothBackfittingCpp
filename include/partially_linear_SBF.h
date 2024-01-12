//
// Created by Asger Morville on 2024/01/12.
// Implements the partially linear smooth backfitting algorithm
//

#ifndef SMOOTH_BACKFITTING_LIBRARY_PARTIALLY_LINEAR_SBF_H
#define SMOOTH_BACKFITTING_LIBRARY_PARTIALLY_LINEAR_SBF_H

#include "smooth_backfitting_core.h"

typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXXd Array;
typedef Eigen::MatrixXd Matrix;

Vector PL_SBF(Vector &Y, Matrix &X, Matrix &Z){
    size_t n = Y.size();
    size_t d_x = X.cols();
    size_t d_z = Z.cols();
    Matrix m_hat_xj_evals(n,d_x);

    AddFunction y_fit_add_func = SBF(Y,Z);
    Vector m_hat_y_evals = y_fit_add_func.predict(Z);

    for (int j = 0; j < d_x; j++){
        Vector colJ = X.col(j);
        AddFunction xjResponse = SBF(colJ, Z);
        m_hat_xj_evals.col(j) = xjResponse.predict(Z);
    }
    Matrix xTilde = X - m_hat_xj_evals;
    Vector yTilde = Y - m_hat_y_evals;
    Matrix A = xTilde.transpose()*xTilde;
    Vector b = xTilde.transpose()*yTilde;
    // Solve Ax = b using Cholesky decomposition,
    Vector betaHat = A.llt().solve(b);
    Vector newY = Y - X*betaHat;
    AddFunction mHatAddFunc = SBF(newY, Z);
    Vector fittedValues = X*betaHat + mHatAddFunc.predict(Z);
    return fittedValues;
}

void plSBFWrapper(double* yPtr, double* xPtr, double* zPtr, double* outputPtr, int n, int d_x, int d_z){
    Eigen::Map<Eigen::VectorXd> yVectorMap(yPtr, n);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> xMatrixMap(xPtr, n, d_x);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> zMatrixMap(zPtr, n, d_z);
    Vector yVector = yVectorMap;
    Matrix xMatrix = xMatrixMap;
    Matrix zMatrix = zMatrixMap;
    Vector fittedValues = PL_SBF(yVector, xMatrix, zMatrix);
    std::copy(fittedValues.data(), fittedValues.data() + fittedValues.size(), outputPtr);
};


#endif //SMOOTH_BACKFITTING_LIBRARY_PARTIALLY_LINEAR_SBF_H

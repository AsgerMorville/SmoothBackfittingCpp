//
// Created by Asger Morville on 2024/01/07.
//

#ifndef SMOOTH_BACKFITTING_LIBRARY_SMOOTH_BACKFITTING_CORE_H
#define SMOOTH_BACKFITTING_LIBRARY_SMOOTH_BACKFITTING_CORE_H

#include <iostream>
#include <Eigen/Dense>
#include "additive_function.h"
#include "quantile.h"
#include <numeric>

typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXXd Array;
typedef Eigen::MatrixXd Matrix;

double epan(double x){
    if (std::abs(x) < 1){
        return (3.0 / 4.0) * (1 - pow(x, 2));
    } else {
        return 0;
    }
}

double kern(double x, double h){
    return (1.0/h)*epan(x/h);
}


double trapz(Vector vector, double dx){
    // This function should
    size_t n = vector.size();
    return (2*vector.sum()-vector[0] - vector[n-1])*dx/2.0;
}

Array createGrid(Vector &xMin, Vector &xMax, int M) {
    size_t d = xMin.size();
    Array grid = Array::Zero(M,d);
    for (int j = 0; j < d; j++) {
        grid.col(j) = Vector::LinSpaced(M, xMin(j), xMax(j));
    }
    return grid;
}

Vector hInitialize(Array &xTranslated){
    std::vector<double> quantileList{0.25,0.75};
    //Vector stdDev = (xTranslated.rowwise()-xTranslated.colwise().mean()).colwise().norm();
    Vector stdDev = sqrt(((xTranslated.rowwise()-xTranslated.colwise().mean()).square()).colwise().mean());
    size_t n = xTranslated.rows();
    size_t d = xTranslated.cols();
    Vector h = Vector::Zero(d);
    for (int j = 0; j < d; j++){
        std::vector<double> columnJ(xTranslated.col(j).data(),xTranslated.col(j).data()+n);
        std::vector<double> quantiles = Quantile(columnJ,quantileList);
        double IQR = quantiles[1] - quantiles[0];
        double h1 = (0.9/pow(n,0.2)) * std::min(stdDev(j),IQR/1.349);

        // Make sure p_j(x) > 0 for all x
        std::sort(columnJ.begin(),columnJ.end());
        std::vector<double> differences(n);
        std::adjacent_difference(columnJ.begin(),columnJ.end(),differences.begin());
        auto maxDifferencePtr = std::max_element(differences.begin()+1,differences.end());
        *maxDifferencePtr *= 0.5;
        h[j] = std::max(h1,*maxDifferencePtr);
    }
    return h;
};

std::vector<Matrix> generateKhTable(Vector &xGrid, Array &xTranslated, Vector &h, double dx){
    size_t M = xGrid.size();
    size_t d = xTranslated.cols();
    size_t n = xTranslated.rows();
    std::vector<Matrix> khTable;
    //khTable.reserve(d);
    for (int j = 0; j < d; j++){
        double hj = h[j];
        // Construct M,n matrix
        Matrix mat1 = Matrix::Zero(M,n);
        for (int i = 0; i < n; i++){
            for (int l=0; l < M; l++){
                mat1(l,i) = kern(xGrid(l)-xTranslated(i,j),hj);
            }
            // Make sure that it integrates to one
            mat1.col(i) /= trapz(mat1.col(i),dx);
        }
        khTable.push_back(mat1);
    }
    return khTable;
}

Array generatePHatTable(std::vector<Matrix> khTable){
    size_t M = khTable[0].rows();
    size_t d = khTable.size();
    Array output = Array::Zero(M,d);
    for (int j=0; j < d; j++){
        for (int l=0; l < M; l++){
            output(l,j) = khTable[j].row(l).mean();
        }
    }
    return output;
}

std::vector<Matrix> generatePHatTable2(std::vector<Matrix> khTable){
    size_t n = khTable[0].cols();
    size_t d = khTable.size();
    std::vector<Matrix> output;
    for (int i=0; i<d; i++){
        for (int j=0; j<d; j++){
            Matrix prod = khTable[i]*khTable[j].transpose()/n;
            output.push_back(prod);
        }
    }
    return output;
}

Array generateFHatTable(Vector yCentered, std::vector<Matrix> khTable, Array pHatTable){
    size_t M = khTable[0].rows();
    size_t n = khTable[0].cols();
    size_t d = khTable.size();
    Array output = Array::Zero(M,d);
    for (int j=0; j<d; j++){
        for (int l=0; l<M; l++){
            output(l,j) = yCentered.dot(khTable[j].row(l))/n;
        }
    }
    output /= pHatTable;
    return output;
}

bool checkConv(Matrix &mHat, Matrix &mOld){
    size_t d = mHat.cols();
    double eps = 0.0001;
    for (int j=0; j<d; j++){
        if ((mHat.col(j)-mOld.col(j)).squaredNorm() / ((mOld.col(j)).squaredNorm()+eps) > eps){
            return false;
        }
    }
    return true;
}

AddFunction SBF(Vector &Y, Array &X){
    // Initialize parameters. M is number of gridpoints.
    size_t M = 100;
    size_t n = X.rows();
    size_t maxIter = 20;
    size_t d = X.cols();

    double yMean = Y.mean();
    Vector yCentered = Y.array()-yMean;
    Vector xGrid = Vector::LinSpaced(M, 0, 1);
    Vector xMinValues = X.colwise().minCoeff();
    Vector xMaxValues = X.colwise().maxCoeff();
    Array xTranslated = Array::Zero(n,d);
    for (int j = 0; j < d; j++){
        xTranslated.col(j) = (X.col(j) - xMinValues(j))/(xMaxValues(j)-xMinValues(j));
    }

    Matrix mHat = Matrix::Zero(M,d);
    double dx = 1.0/(M-1.0);

    Vector h = hInitialize(xTranslated);
    std::vector<Matrix> khTable = generateKhTable(xGrid, xTranslated, h, dx);
    Array pHatTable = generatePHatTable(khTable);
    std::vector<Matrix> pHatTable2 = generatePHatTable2(khTable);
    Array fHatTable = generateFHatTable(yCentered,khTable,pHatTable);

    // Optimization loop
    for (int B = 0; B < maxIter; B++){
        Matrix mOld = mHat;
        for (int j = 0; j<d; j++){
            for (int l = 0; l < M; l++){
                double integral_sum = 0;
                for (int k = 0; k < d; k++){
                    if (k != j){
                        double acc = mHat.col(k).dot(pHatTable2[j*d+k].row(l));
                        integral_sum += 2*acc - mHat(0,k)*pHatTable2[d*j+k](l,0)-mHat(M-1,k)*pHatTable2[j*d+k](l,M-1);
                    }
                }
                mHat(l, j) = fHatTable(l, j) - (integral_sum * dx) / (2*pHatTable(l, j));
            }
        }
        if (checkConv(mHat, mOld)){
            break;
        }
    }
    Array fin_x_grid = createGrid(xMinValues, xMaxValues, M);
    AddFunction output = AddFunction(fin_x_grid,mHat, yMean);
    return output;
};


void sbfWrapper(double* yPtr, double* xPtr, double* outputPtr, int n, int d){
    Eigen::Map<Eigen::VectorXd> yVectorMap(yPtr, n);
    Eigen::Map<Eigen::ArrayXXd> xArrayMap(xPtr, n, d);
    Vector yVector = yVectorMap;
    Array xArray = xArrayMap;
    AddFunction output = SBF(yVector, xArray);
    Vector fittedValues = output.predict(xArray);
    //outputPtr = fittedValues.data();
    std::copy(fittedValues.data(), fittedValues.data() + fittedValues.size(), outputPtr);
};


/*

class sbfFitter{
public:
    sbfFitter () = default;
    AddFunction fit(Vector &Y, Array &X){
        return SBF(Y,X);
    }
};




Eigen::ArrayXXd SBF_new(Eigen::VectorXd Y, Eigen::ArrayXXd X,const size_t n)
{
    const size_t M = 100;
    size_t maxIter =20;
    size_t d = X.cols();
    //std::cout << "d" << d << "\n";
    //const size_t n = X.rows();

    double Y_mean = Y.mean();
    Eigen::VectorXd Y_center = Y.array()- Y_mean;
    Eigen::VectorXd x_grid = Eigen::VectorXd::LinSpaced(M,0,1);
    Eigen::ArrayXd X_min = X.colwise().minCoeff();
    Eigen::ArrayXd X_max = X.colwise().maxCoeff();
    Eigen::ArrayXXd X_center = Eigen::ArrayXXd::Zero(n,d);
    for (int j = 0; j < d; j++){
        X_center.col(j) = (X.col(j) - X_min(j))/(X_max(j)-X_min(j));
    }

    Eigen::MatrixXd m_hat = Eigen::MatrixXd::Zero(M,d);
    double dx = 1.0/(M-1);

    auto t00 = std::chrono::steady_clock::now();
    //h initialization
    std::vector<double> h(d,0);
    std::vector<double> quantile_list{0.25,0.75};
    auto std_dev = sqrt(((X_center.rowwise()-X_center.colwise().mean()).square()).colwise().mean());
    //std::cout << "STD DEV" << std_dev << "\n";
    for (int j = 0; j < d; j++){
        std::vector<double> xx(X_center.col(j).data(),X_center.col(j).data()+n);
        std::vector<double> qq = Quantile(xx,quantile_list);
        double IQR = qq[1] - qq[0];
        double h1 = (0.9/pow(n,0.2)) * std::min(std_dev(0,j),IQR/1.349);

        // Make sure p_j(x) > 0 for all x
        std::sort(xx.begin(),xx.end());
        std::vector<double> differences(n);
        std::adjacent_difference(xx.begin(),xx.end(),differences.begin());
        auto maxim_elem = std::max_element(differences.begin()+1,differences.end());
        *maxim_elem *= 0.5;
        double min_h = std::max(*maxim_elem,h[j])+0.001;
        h[j] = std::max(h1,min_h);
    }

    auto t01 = std::chrono::steady_clock::now();
    //K_table initialization
    std::vector<Eigen::MatrixXd> Kh_table;
    Kh_table.reserve(d);
    for (int j = 0; j < d; j++){
        double hj = h[j];
        // Construct M,n matrix
        Eigen::ArrayXXd mat1 = Eigen::ArrayXXd::Zero(M,n);
        for (int i = 0; i < n; i++){
            for (int l=0; l < M; l++){
                mat1(l,i) = kern(x_grid(l)-X_center(i,j),hj);
            }
            // Make sure that it integrates to one
            mat1.col(i) /= trapz1(mat1.col(i),dx);
            //std::cout << "Norm constant: " <<normalizing_const << "\n";
        }
        Kh_table.push_back(mat1);
    }

    //p_hat_table initialization
    auto t02 = std::chrono::steady_clock::now();
    Eigen::ArrayXXd p_hat_table = Eigen::ArrayXXd::Zero(M,d);
    for (int j=0; j < d; j++){
        for (int l = 0; l< M; l++){
            //std::cout << "Kh_table[j].row(l)" << Kh_table[j].row(l) << "\n";
            p_hat_table(l,j) = Kh_table[j].row(l).mean();
        }
    }
    //std::cout << "check" << "\n";
    //p_hat_table2 initialization
    auto t03 = std::chrono::steady_clock::now();
    std::vector<Eigen::MatrixXd> p_hat_table2;
    p_hat_table2.reserve(d*d);
    for (int i = 0; i < d; i++){
        for (int j = 0; j < d; j++){
            Eigen::MatrixXd prod = Kh_table[i]*Kh_table[j].transpose()/n;
            p_hat_table2.push_back(prod);
        }
    }
    //std::cout << "check" << "\n";
    //f_hat_table initialization
    auto t04 = std::chrono::steady_clock::now();
    Eigen::ArrayXXd f_hat_table = Eigen::ArrayXXd::Zero(M,d);
    for (int j=0; j<d; j++){
        for (int l=0; l<M; l++){
            f_hat_table(l,j) = Y_center.dot(Kh_table[j].row(l))/n;
        }
    }
    f_hat_table /= p_hat_table;

    auto t05 = std::chrono::steady_clock::now();
    //Optimization loop.

    for (int B = 0; B < maxIter; B++){
        for (int j = 0; j<d; j++){
            for (int l = 0; l < M; l++){
                double integral_sum = 0;
                for (int k = 0; k < d; k++){
                    if (k != j){
                        double acc = m_hat.col(k).dot(p_hat_table2[j*d+k].row(l));
                        integral_sum += 2*acc - m_hat(0,k)*p_hat_table2[d*j+k](l,0)-m_hat(M-1,k)*p_hat_table2[j*d+k](l,M-1);
                        //integral_sum += 2*acc;
                    }
                }
                m_hat(l, j) = f_hat_table(l, j) - (integral_sum * dx) / (2*p_hat_table(l, j));
            }
        }
    }

    auto t06 = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::milli> h_init_time = t01-t00;
    std::chrono::duration<double, std::milli> Kh_init_time = t02-t01;
    std::chrono::duration<double, std::milli> phat_init_time = t03-t02;
    std::chrono::duration<double, std::milli> phat2_init_time = t04-t03;
    std::chrono::duration<double, std::milli> f_hat_init_time = t05-t04;
    std::chrono::duration<double, std::milli> opt_loop_time = t06-t05;

    std::cout << "h init time: " << h_init_time.count() << " ms" << "\n";
    std::cout << "Kh init time: " << Kh_init_time.count()<< " ms" << "\n";
    std::cout << "phat init time: " << phat_init_time.count()<< " ms"  << "\n";
    std::cout << "phat2 init time: " << phat2_init_time.count()<< " ms" << "\n";
    std::cout << "fhat init time: " << f_hat_init_time.count()<< " ms"<< "\n";
    std::cout << "opt loop time: " << opt_loop_time.count()<< " ms"  << "\n";
    return m_hat;
}
*/
#endif //SMOOTH_BACKFITTING_LIBRARY_SMOOTH_BACKFITTING_CORE_H

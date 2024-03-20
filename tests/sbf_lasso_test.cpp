//
// Created by Asger Morville on 2024/02/23.
//
#import "smooth_backfitting_core.h"

int main(){
    int n = 100;
    int d = 10;
    Matrix X = Matrix::Random(n,d);
    Vector Y = X.col(0).array().square();
    std::cout << Y << "\n";
    AddFunction add_func_test = SBFLasso(Y,X,0.01);

    return 0;
};
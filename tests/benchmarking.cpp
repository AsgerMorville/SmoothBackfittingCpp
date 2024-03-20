//
// Created by Asger Morville on 2024/01/28.
//

#import "smooth_backfitting_core.h"

int main(){
    int n = 2000;
    int d = 20;
    Matrix X = Matrix::Random(n,d);
    Vector Y = Vector::Random(n);

    AddFunction add_func_test = SBF(Y,X);

    return 0;
};
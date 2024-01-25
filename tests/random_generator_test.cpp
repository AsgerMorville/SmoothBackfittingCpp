//
// Created by Asger Morville on 2024/01/25.
//
#include "utils/random_generator.h"
#include <iostream>
int main(){
    int n = 10;
    int d = 2;

    std::mt19937 mt(1);

    Matrix Test = uniformMatrix(n,d,mt);
    Vector Test2 = stdNormalVector(n, mt);
    Matrix bernoulliTest = bernoulliMatrix(n,d,mt);
    std::cout << bernoulliTest << "\n";

    return 0;
}
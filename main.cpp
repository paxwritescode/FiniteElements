#include <iostream>
#include "GalerkinMethod.h"

int main(void){
    std::cout << "Hello lab1" << std::endl;
    double* c = new double[4] {1, 2, -4, -8};
    double* b = new double[5] {2, 9, 17, 15, 3};
    double* a = new double[4] {2, 4, 4, 2};
    double* d = new double[5] {-10, -26, -16, -2, 16};

    double* x = new double[5];

    TridiagonalMatrixAlgorithm(5, a, b, c, d, x);

    for (int i = 0; i < 5; i++)
    {
        std::cout << x[i] << std::endl;
    }
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
}

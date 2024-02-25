#include <cmath>
#include <iostream>
#include "GalerkinMethod.h"
#include "RitzMethod.h"

double TestFunc(double x)
{
    return sin(x);
}

void TestSimpson()
{
    printf("%lf - Simpson computed\n", SimpsonIntegrate(100, TestFunc, PhiParams{0, PI/2, PI}, 0, PI));
    printf("%lf - exact\n\n", 4/PI);
}

void TestTridiagonal()
{
    double* b = new double[5] {0, 2, 4, 4, 2};
    double* c = new double[5] {2, 9, 17, 15, 3};
    double* d = new double[5] {1, 2, -4, -8, 0};
    double* r = new double[5] {-10, -26, -16, -2, 16};

    double* x = new double[5];

    TridiagonalMatrixAlgorithm(5, r, b, c, d, x);

    for (int i = 0; i < 5; i++)
    {
        printf("%lf ", x[i]);
    }

    printf("\nCorrect:\n");
    printf("-4, -2, 0, 2, 4\n\n");
    
    delete[] b;
    delete[] c;
    delete[] d; 
    delete[] r;   
}

void TestRitzFunctional()
{
    std::cout << "TestRitzFunctional" << std::endl;
    int N = 5;
    double a = -1, b = 1;
    double* c = new double[N + 1];
    for (int i = 0; i < N + 1; i++)
        c[i] = 1;

    double J1 = ComputeJ1(N, c, a, b);

    std::cout << "J1 - calculated value: " << J1 << std::endl;
    std::cout << "J1 - correct value: " << 0 << std::endl;

    double J2 = ComputeJ2(N, c, f, a, b);

    std::cout << "J2 - calculated value: " << J2 << std::endl;
    std::cout << "J2 - correct value: " << 2 << std::endl;

    delete[] c;
}

void TestMinBetween3Numbers()
{
    double a = .5, b = 1.5, c = -.5;
    std::cout << MinBetween3Numbers(a, b, c) << std::endl;
    std::cout << MinBetween3Numbers(b, c, a) << std::endl;
    std::cout << MinBetween3Numbers(c, b, a) << std::endl;
    std::cout << MinBetween3Numbers(c, a, b) << std::endl;
}


void Test()
{
    TestSimpson();
    TestTridiagonal();
    TestRitzFunctional();
    TestMinBetween3Numbers();
}

int main(void)
{
    Test();
}
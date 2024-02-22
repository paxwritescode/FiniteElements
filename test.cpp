#include "test.h"

double TestFunc(double x)
{
    return sin(x);
}

void TestSimpson()
{
    printf("%lf - Simpson computed\n", SimpsonIntegrate(100, TestFunc, PhiParams{0, PI/2, PI}, 0, PI));
    printf("%lf - exact\n", 4/PI);
}

void TestTridiagonal()
{
    double* b = new double[10] {0, -5, -5, -5, -5, -5, -5, -5, -5, -5};
    double* c = new double[10] {5, 10, 10, 10, 10, 10, 10, 10, 10, 5};
    double* d = new double[10] {-5, -5, -5, -5, -5, -5, -5, -5, -5, 0};
    double* r = new double[10] {0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1};

    double* x = new double[10];

    TridiagonalMatrixAlgorithm(10, r, b, c, d, x);

    for (int i = 0; i < 10; i++)
    {
        printf("%lf ", x[i]);
    }

    // printf("\nCorrect:\n");
    // printf("-4, -2, 0, 2, 4");

    printf("\n\n");
    
    delete[] b;
    delete[] c;
    delete[] d; 
    delete[] r;   
}

void Test()
{
    TestSimpson();
    TestTridiagonal();
}
#include "test.h"

double TestFunc(double x)
{
    return sin(x);
}

void TestSimpson()
{
    printf("%lf - Simpson computed\n", SimpsonIntegrate(100, TestFunc, PhiParams{0, PI/2, PI}));
    printf("%lf - exact\n", 4/PI);
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
    printf("-4, -2, 0, 2, 4");
    
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
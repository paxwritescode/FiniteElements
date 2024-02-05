#include <iostream>
#include <iomanip>

void TridiagonalMatrixAlgorithm(int n, double *a, double *b, double *c, double *d, double *x)
{
    //a - lower diagonal
    //b - main diagonal
    //c - upper diagonal
    //d - right hand side

    // Forward elimination
    for (int i = 1; i < n; ++i)
    {
        double factor = a[i] / b[i - 1];
        b[i] -= factor * c[i - 1];
        d[i] -= factor * d[i - 1];
    }

    // Backward substitution
    x[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
}

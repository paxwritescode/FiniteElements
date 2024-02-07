#include <cstdio>
#include <cstdlib>

void TridiagonalMatrixAlgorithm(int n, double *r, double *b, double *c, double *d, double *x)
{
    double* delta = (double*)calloc(n, sizeof(double));
    double* lambda = (double*)calloc(n, sizeof(double));

    delta[0] = -d[0] / c[0];
    lambda[0] = r[0] / c[0];

    for (int i = 1; i < n; i++)
    {
        delta[i] = -d[i] / (c[i] + b[i] * delta[i - 1]);
        lambda[i] = (r[i] - b[i] * lambda[i - 1]) / (c[i] + delta[i - 1] * b[i]);
    }

    x[n - 1] = lambda[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = delta[i] * x[i + 1] + lambda[i];
    }

    free(delta);
    free(lambda);
}

#include "GalerkinMethod.h"
#include "functions.h"

void TridiagonalMatrixAlgorithm(int n, double *r, double *b, double *c, double *d, double *x)
{
    double *delta = (double *)calloc(n, sizeof(double));
    double *lambda = (double *)calloc(n, sizeof(double));

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

void GalerkinMethod(double a, double b, int n, double *uj)
{
    // Constructing tridiagonal Matrix
    double *d = (double *)calloc(n, sizeof(double));  // a_(j - 1)(j), upper diagonal
    double *c = (double *)calloc(n, sizeof(double));  // a_(j)(j), main diagonal
    double *bm = (double *)calloc(n, sizeof(double)); // a(j + 1)(j), lower diagonal
    double *r = (double *)calloc(n, sizeof(double));

    bm[0] = 0, d[n - 1] = 0;

    for (int j = 1; j < n - 1; j++)
    {
        PhiParams phiParams = FillPhiParams(a, b, n, j);

        c[j] = 1 / (phiParams.x_i - phiParams.x_p) + 1 / (phiParams.x_n - phiParams.x_i);
        if (j != 0)
            bm[j] = -1 / (phiParams.x_n - phiParams.x_i);
        if (j != n - 1)
            d[j] = -1 / (phiParams.x_i - phiParams.x_p);

        r[j] = SimpsonIntegrate(100, f, phiParams, a, b);
    }

    c[0] = 1 / (ComputeRegularGridNode(a, b, 1, n) - a);
    c[n - 1] = 1 / (b - ComputeRegularGridNode(a, b, n - 1, n));
    r[0] = SimpsonIntegrate0(20, f, {0, a, ComputeRegularGridNode(a, b, 1, n)}, a, b);
    //r[n - 1] = SimpsonIntegraten(20, f, {ComputeRegularGridNode(a, b, n - 2, n), b, 0}, a, b);
    r[n - 1] = 0.1;
    bm[n - 1] = -1 / (b - ComputeRegularGridNode(a, b, n - 1, n));
    d[0] = -1 / (ComputeRegularGridNode(a, b, 1, n) - a);

    TridiagonalMatrixAlgorithm(n, r, bm, c, d, uj);

    free(d);
    free(c);
    free(bm);
    free(r);

    //return uj;
}


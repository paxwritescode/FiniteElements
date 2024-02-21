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

    for (int j = 0; j < n; j++)
    {
        PhiParams phiParams = FillPhiParams(a, b, n, j);

        // TODO boundaries

        c[j] = 1 / (phiParams.x_i - phiParams.x_p) + 1 / (phiParams.x_n - phiParams.x_i);
        if (j != 0)
            bm[j] = -1 / (phiParams.x_n - phiParams.x_i);
        if (j != n - 1)
            d[j] = -1 / (phiParams.x_i - phiParams.x_p);

        r[j] = SimpsonIntegrate(100, f, phiParams, a, b);
    }

    TridiagonalMatrixAlgorithm(n, r, bm, c, d, uj);

    free(d);
    free(c);
    free(bm);
    free(r);

    //return uj;
}
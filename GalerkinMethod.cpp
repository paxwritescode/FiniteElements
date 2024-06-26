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
    double *d = (double *)calloc(n + 1, sizeof(double));  // a_(j - 1)(j), upper diagonal
    double *c = (double *)calloc(n + 1, sizeof(double));  // a_(j)(j), main diagonal
    double *bm = (double *)calloc(n + 1, sizeof(double)); // a(j + 1)(j), lower diagonal
    double *r = (double *)calloc(n + 1, sizeof(double));

    bm[0] = 0, d[n] = 0;

    for (int j = 1; j < n; j++)
    {
        PhiParams phiParams = FillPhiParams(a, b, n, j);

        c[j] = 1 / (phiParams.x_i - phiParams.x_p) + 1 / (phiParams.x_n - phiParams.x_i);
        if (j != 0)
            bm[j] = -1 / (phiParams.x_n - phiParams.x_i);
        if (j != n)
            d[j] = -1 / (phiParams.x_i - phiParams.x_p);

        r[j] = SimpsonIntegrate(200, f, phiParams, a, b);
    }

    c[0] = 1 / (ComputeRegularGridNode(a, b, 1, n) - a);
    c[n] = 1 / (b - ComputeRegularGridNode(a, b, n - 1, n));

    // printf("\nmain diagonal:\n");
    // for (int i = 0; i <= n; i++)
    //     printf("%lf ", c[i]);
    // printf("\n");

    // printf("\nupper diagonal d:\n");
    // for (int i = 0; i <= n; i++)
    //     printf("%lf ", d[i]);
    // printf("\n");   

    // printf("\nlower diagonal b:\n");
    // for (int i = 0; i <= n; i++)
    //     printf("%lf ", bm[i]);
    // printf("\n");  

    // printf("\nright column r:\n");
    // for (int i = 0; i <= n; i++)
    //     printf("%lf ", r[i]);
    // printf("\n");

    TridiagonalMatrixAlgorithm(n + 1, r, bm, c, d, uj);

    // printf("\narray of coefficients u_j\n");
    // for (int i = 0; i <= n; i++)
    //     printf("%lf ", uj[i]);
    // printf("\n\n");

    free(d);
    free(c);
    free(bm);
    free(r);

    //return uj;
}
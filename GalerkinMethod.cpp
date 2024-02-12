#include "GalerkinMethod.h"

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

double ComputeRegularGridNode(double a, double b, int i, int n)
{
    double h = (b - a) / n;
    return a + i * h;
}

double ComputeAdaptiveGridNode(double a, double b, int i, int n)
{
    // TODO choose a function
    double h = (b - a) / n;
    return a + i * h;
}

double ComputeBasicFunction(int i, double x, double x_p, double x_i, double x_n) // x_p = x_(i - 1), previous; x_n = x_(i + 1), next
{
    if (x_p < x && x_i >= x)
        return (x - x_p) / (x_i - x_p);
    else if (x_i < x && x_n >= x)
        return (x_n - x) / (x_n - x_i);
    else
        return 0;
}


double phi(double x, PhiParams params){
    return ComputeBasicFunction(params.i, x, params.x_p, params.x_i, params.x_n);
}


double ComputeTestFunction(int i, double x, double x_p, double x_i, double x_n) // x_p = x_(i - 1), previous; x_n = x_(i + 1), next
{
    if (x_p < x && x_i >= x)
        return 1 / (x_i - x_p);
    else if (x_i < x && x_n >= x)
        return -1 / (x_n - x_i);
    else
        return 0;
}

double SimpsonIntegrate(int n, double (*func)(double), PhiParams phiParams) //TODO 76 line f * phi
{
    const double length = (phiParams.x_n - phiParams.x_p) / n;

    double simpson_integral = 0;
    for (int step = 0; step < n; step++)
    {
        const double x1 = phiParams.x_p + step * length;
        const double x2 = phiParams.x_p + (step + 1) * length;
        simpson_integral += (x2 - x1) / 6.0 * (func(x1) * phi(x1, phiParams) + 
        4.0 * func(0.5 * (x1 + x2)) * phi(0.5 * (x1 + x2), phiParams) + func(x2) * phi(x2, phiParams));
    }

    return simpson_integral;
}

double GalerkinMethod(double x, double a, double b, int n)
{
    double u_numeric = 0;

    // Constructing tridiagonal Matrix
    double *d = (double *)calloc(n, sizeof(double));  // a_(j - 1)(j), upper diagonal
    double *c = (double *)calloc(n, sizeof(double));  // a_(j)(j), main diagonal
    double *bm = (double *)calloc(n, sizeof(double)); // a(j + 1)(j), lower diagonal
    double *uj = (double *)calloc(n, sizeof(double));
    double *r = (double *)calloc(n, sizeof(double));

    bm[0] = 0, d[n - 1] = 0;

    for (int j = 0; j < n; j++)
    {
        PhiParams phi = {0};
        phi.x_i = ComputeRegularGridNode(a, b, j, n);
        phi.x_p = ComputeRegularGridNode(a, b, j - 1, n); // previous
        phi.x_n = ComputeRegularGridNode(a, b, j + 1, n); // next
        phi.i = j;

        c[j] = 1 / ((phi.x_i - phi.x_p) * (phi.x_i - phi.x_p)) 
        + 1 / ((phi.x_n - phi.x_i) * (phi.x_n - phi.x_i));
        if (j != 0)
            bm[j] = -(1 / (phi.x_n - phi.x_i) * (phi.x_n - phi.x_i));
        if (j != n - 1)
            d[j] = -(1 / (phi.x_i - phi.x_p) * (phi.x_i - phi.x_p));

        r[j] = SimpsonIntegrate(20, f, phi);
    }


    TridiagonalMatrixAlgorithm(n, r, bm, c, d, uj);

    for (int j = 0; j < n; j++)
    {
        PhiParams phi = {0};
        phi.x_i = ComputeRegularGridNode(a, b, j, n);
        phi.x_p = ComputeRegularGridNode(a, b, j - 1, n); // previous
        phi.x_n = ComputeRegularGridNode(a, b, j + 1, n); // next
        phi.i = j;

        u_numeric += uj[j] * ComputeBasicFunction(j, x, phi.x_p, phi.x_i, phi.x_n); //*phi_j(x)
    }

    free(d);
    free(c);
    free(bm);
    free(uj);
    free(r);

    return u_numeric;
}
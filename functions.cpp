#include <cmath>
#include "functions.h"

double f(double x)
{
    return - PI * cos(PI * x / 2) + PI * PI / 4 * x * sin(PI * x / 2);
    //return 1;
}

double u_exact(double x)
{
    return x * sin(PI * x / 2);
    //return -x * x / 2 + 0.5;
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

double phi(double x, PhiParams params, double a, double b)
{
    // x_p = x_(i - 1), previous; x_n = x_(i + 1), next
    if (x < a || x > b)
        return 0;
    if (params.x_p < x && params.x_i >= x)
        return (x - params.x_p) / (params.x_i - params.x_p);
    else if (params.x_i < x && params.x_n >= x)
        return (params.x_n - x) / (params.x_n - params.x_i);
    else
        return 0;
}

double der_phi(double x, PhiParams params, double a, double b)
{
    // x_p = x_(i - 1), previous; x_n = x_(i + 1), next
    if (x < a || x > b)
        return 0;
    if (params.x_p < x && params.x_i >= x)
        return 1 / (params.x_i - params.x_p);
    else if (params.x_i < x && params.x_n >= x)
        return -1 / (params.x_n - params.x_i);
    else
        return 0;
}

double SimpsonIntegrate(int n, double (*func)(double), PhiParams phiParams, double a, double b)
{
    const double length = (phiParams.x_n - phiParams.x_p) / n;

    double simpson_integral = 0;
    for (int step = 0; step < n; step++)
    {
        const double x1 = phiParams.x_p + step * length;
        const double x2 = phiParams.x_p + (step + 1) * length;
        simpson_integral += (x2 - x1) / 6.0 * (func(x1) * phi(x1, phiParams, a, b) + 4.0 * func(0.5 * (x1 + x2)) * phi(0.5 * (x1 + x2), phiParams, a, b) + func(x2) * phi(x2, phiParams, a, b));
    }

    return simpson_integral;
}

PhiParams FillPhiParams(double a, double b, int n, int j)
{
    PhiParams phiParams = {0};

    phiParams.x_i = ComputeRegularGridNode(a, b, j, n);
    phiParams.x_p = ComputeRegularGridNode(a, b, j - 1, n); // previous
    phiParams.x_n = ComputeRegularGridNode(a, b, j + 1, n); // next

    return phiParams;
}

double ComputeErrorNorm(double* array1, double* array2, int n)
{
    double norm = 0;
    for(int i = 0; i < n; i++)
    {
        norm += sqrt((array1[i] - array2[i]) * (array1[i] - array2[i]));
    }

    return norm;
}

double ComputeSolution(double x, double a, double b, int n, double *uj)
{
    double u_numeric = 0;

    for (int j = 0; j < n; j++)
    {
        PhiParams phiParams = FillPhiParams(a, b, n, j);

        u_numeric += uj[j] * phi(x, phiParams, a, b);
    }

    return u_numeric;
}
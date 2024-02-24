#include "RitzMethod.h"
#include "functions.h"

/** Calculate integral from a to b of (u')^2
 * @param N: number of basic functions
 * @param c: array of coefficients c_i, i from 0 to N-1
 * @param a: left boundary
 * @param b: right boundary
 */
double ComputeJ1(const int N, const double *const c, const double a, const double b)
{
    double res = 0;
    // main diagonal
    for (int j = 1; j < N - 1; j++)
    {
        double x_p = ComputeRegularGridNode(a, b, j - 1, N); // x_{j - 1}
        double x_j = ComputeRegularGridNode(a, b, j, N);
        double x_n = ComputeRegularGridNode(a, b, j + 1, N); // x_{j + 1}
        res += c[j] * c[j] * (1 / (x_j - x_p) + 1 / (x_n - x_j));
    }
    res += c[0] * c[0] * (1 / ComputeRegularGridNode(a, b, 1, N) - a);
    res += c[N - 1] * c[N - 1] * (b - ComputeRegularGridNode(a, b, N - 1, N));

    // upper diagonal
    for (int j = 1; j < N; j++)
    {
        double x_p = ComputeRegularGridNode(a, b, j - 1, N); // x_{j - 1}
        double x_j = ComputeRegularGridNode(a, b, j, N);
        res += c[j - 1] * c[j] * (-1 / (x_j - x_p));
    }

    // lower diagonal
    for (int j = 0; j < N - 1; j++)
    {
        double x_j = ComputeRegularGridNode(a, b, j, N);
        double x_n = ComputeRegularGridNode(a, b, j + 1, N); // x_{j + 1}
        res += c[j] * c[j + 1] * (-1 / (x_n - x_j));
    }

    return res;
}

/** Calculate integral from a to b of {f(x) * u(x)}
 * computing as sum_{j = 0}^{N - 1} c_j \int_{x_j}^{x_{j + 1}} f(x) phi_j(x) dx
 * @param N: number of basic functions
 * @param c: array of coefficients c_i, i from 0 to N-1
 * @param func: f (function of right part of equaiion)
 * @param a: left boundary
 * @param b: right boundary
 *
 */
double ComputeJ2(const int N, const double *const c, double(*func)(double), const double a, const double b)
{
    double res = 0;
    for (int j = 0; j < N; j++)
    {
        PhiParams phiParams = FillPhiParams(a, b, N, j);
        res += c[j] * SimpsonIntegrate(20, func, phiParams, phiParams.x_i, phiParams.x_n);
    }

    return res;
}

double ComputeFunctionalValue(int N, double *c, double (*func)(double), double a, double b)
{
    // J1 = int_a^b u'^2 dx
    double J1 = ComputeJ1(N, c, a, b);
    // J2 = int_a^b f * u dx
    double J2 = ComputeJ2(N, c, func, a, b);

    double J = J1 - J2;

    return J;
}

double MinBetween3Numbers(double a, double b, double c)
{
    double min = a;
    if (min > b)
        min = b;
    if (min > c)
        min = c;
    return min;
}

double *RitzMethod(int N, double a, double b) // N is the number of basic functions
{
    // RETURNS: an array of coefficients c_j
    double *c = (double *)calloc(N, sizeof(double));
    // init the array of coefficients by random values from -1 to 1
    for (int j = 0; j < N; j++)
    {
        int random_integer = rand();
        double random_value = ((double)random_integer / RAND_MAX) * 2.0 - 1.0;
        c[j] = random_value;
    }

    double *c_prev = (double *)calloc(N, sizeof(double));
    int CountOfIterations = 0;
    int CountOfIterations_Max = 1e6;

    double* deltas = new double[N];
    for (int i = 0; i < N; i++)
        deltas[i] = DELTA_START;

    while (ComputeErrorNorm(c, c_prev, N) > EPS && CountOfIterations <= CountOfIterations_Max)
    {
        for (int j = 0; j < N; j++)
            c_prev[j] = c[j];

        int j = rand() % N;

        // 1) c_j - delta
        c[j] -= deltas[j];
        double J1 = ComputeFunctionalValue(N, c, f, a, b);
        // 2) c_j
        c[j] += deltas[j];
        double J2 = ComputeFunctionalValue(N, c, f, a, b);
        // 3) c_j + delta
        c[j] += deltas[j];
        double J3 = ComputeFunctionalValue(N, c, f, a, b);

        if (MinBetween3Numbers(J1, J2, J3) == J1)
            c[j] -= 2 * deltas[j];
        else if (MinBetween3Numbers(J1, J2, J3) == J2) // we need to decrease delta
        {
            c[j] -= deltas[j];
            deltas[j] /= 2;
        }
        else
            continue;

        CountOfIterations++;
    }

    if (CountOfIterations > CountOfIterations_Max)
        printf("The method did not converge!!!\n\n");

    free(c_prev);
    delete[] deltas;

    return c;
}
#include "RitzMethod.h"
#include "functions.h"

double ComputeFunctionalValue(int N, double *c, double (*func)(double))
{
    return 0;
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
    //RETURNS: an array of coefficients c_j
    double *c = (double *)calloc(N, sizeof(double));
    // init the array of coefficients by random values from -1 to 1
    for (int j = 0; j < N; j++)
    {
        int random_integer = rand();
        double random_value = ((double)random_integer / RAND_MAX) * 2.0 - 1.0;
        c[j] = random_value;
    }
    double *c_prev = (double *)calloc(N, sizeof(double));
    while (ComputeErrorNorm(c, c_prev, N) > EPS)
    {
        for (int j = 0; j < N; j++)
            c_prev[j] = c[j];

        int j = rand() % N;

        // 1) c_j - delta
        c[j] -= DELTA;
        double J1 = ComputeFunctionalValue(N, c, f);
        // 2) c_j
        c[j] += DELTA;
        double J2 = ComputeFunctionalValue(N, c, f);
        // 3) c_j + delta
        c[j] += DELTA;
        double J3 = ComputeFunctionalValue(N, c, f);

        if (MinBetween3Numbers(J1, J2, J3) == J1)
            c[j] -= 2 * DELTA;
        else if (MinBetween3Numbers(J1, J2, J3) == J2)
            c[j] -= DELTA;
        else
            continue;
    }
    free(c_prev);
    return c;
}
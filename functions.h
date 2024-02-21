#pragma once

#define _A_ 0
#define _B_ 2
#define PI 3.1415926535897932
#define EPS 1e-3
#define DELTA 1e-4

struct PhiParams
{
    double x_p;
    double x_i;
    double x_n;
};
typedef struct PhiParams PhiParams;
 
double f(double x);
double u_exact(double x);
double ComputeRegularGridNode(double a, double b, int i, int n);
double ComputeAdaptiveGridNode(double a, double b, int i, int n);
double phi(double x, PhiParams params, double a, double b);
double der_phi(double x, PhiParams params, double a, double b);
double SimpsonIntegrate(int n, double (*func)(double), PhiParams phiParams, double a, double b);
PhiParams FillPhiParams(double a, double b, int n, int j);
double ComputeErrorNorm(double* array1, double* array2, int n);
double ComputeSolution(double x, double a, double b, int n, double *uj);
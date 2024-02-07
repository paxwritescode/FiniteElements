#pragma once

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "functions.h"

struct PhiParams{
    int i;
    double x_p;
    double x_i;
    double x_n;
};
typedef struct PhiParams PhiParams;

void TridiagonalMatrixAlgorithm(int n, double *r, double *b, double *c, double *d, double *x);
double SimpsonIntegrate(double a, double b, int n, PhiParams phiParams);
double GalerkinMethod(double x, double a, double b, int n);
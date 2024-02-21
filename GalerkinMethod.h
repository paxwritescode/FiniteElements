#pragma once

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "functions.h"

void TridiagonalMatrixAlgorithm(int n, double *r, double *b, double *c, double *d, double *x);
double SimpsonIntegrate(int n, double (*func)(double), PhiParams phiParams, double a, double b);
void GalerkinMethod(double a, double b, int n, double *uj);
double ComputeSolution(double x, double a, double b, int n, double *uj);
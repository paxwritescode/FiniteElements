#pragma once

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "functions.h"

void TridiagonalMatrixAlgorithm(int n, double *r, double *b, double *c, double *d, double *x);
double SimpsonIntegrate(double a, double b, double (*f)(double), int n);
#include <stdlib.h>
#include <cmath>

double ComputeJ1(const int N, const double *const c, const double a, const double b);
double ComputeJ2(const int N, const double *const c, double(*func)(double), const double a, const double b);
double ComputeFunctionalValue(int N, double *c, double (*func)(double), double a, double b);
double MinBetween3Numbers(double a, double b, double c);
double *RitzMethod(int N, double a, double b);
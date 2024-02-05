#include <cmath>
#include "functions.h"

#define PI 3.1415926535897932

double f(double x)
{
    return PI * cos(PI * x / 2) - PI * PI / 4 * x * sin(PI * x / 2);
}

double u_exact(double x)
{
    return x * sin(PI * x / 2);
}
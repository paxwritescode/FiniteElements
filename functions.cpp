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
#include "test.h"

double TestFunc(double x)
{
    return sin(x);
}

void TestSimpson()
{
    printf("%lf - Simpson computed\n", SimpsonIntegrate(100, TestFunc, PhiParams{0, 0, PI/2, PI}));
    printf("%lf - exact\n", 4/PI);
}
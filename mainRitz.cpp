#include <iostream>
#include "RitzMethod.h"
#include "functions.h"

int main(void)
{
    int N = 10;
    printf("Ritz method, %d nodes\n\n", N);
    printf("%d nodes\n\n", N + 1);
    double *c = RitzMethod(N, _A_, _B_);
    for (double x = _A_; x <= _B_ + 0.0001; x += 0.1)
    {
        printf("%lf, ", ComputeSolution(x, _A_, _B_, N, c));
    }
    
    free(c);
    return 0;
}
#include <iostream>
#include "GalerkinMethod.h"
#include "functions.h"
#include "test.h"

int main(void)
{
    int n;
    printf("Input a number of nodes: ");
    scanf("%d", &n);

    double* uj = (double*)calloc(n, sizeof(double));
    GalerkinMethod(_A_, _B_, n, uj);

    for (double x = _A_; x <= _B_ + 0.0001; x += 0.1)
    {
        printf("%lf, ", ComputeGalerkinSolution(x, _A_, _B_, n, uj));
    }

    printf("\n\n");

    for (double x = _A_; x <= _B_ + 0.0001; x += 0.1)
    {
        printf("%lf, ", x);
    }

    // printf("\n");

    // Test();
    
    free(uj);
    return 0;
}

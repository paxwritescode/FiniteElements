#include <iostream>
#include "GalerkinMethod.h"
#include "functions.h"

int main(void)
{
    int n = 10;
    // printf("Input a number of nodes: ");
    // scanf("%d", &n);
    printf("Galerkin method, %d nodes\n\n", n);

    double* uj = (double*)calloc(n + 1, sizeof(double));
    GalerkinMethod(_A_, _B_, n, uj);

    printf("y = ");
    for (double x = _A_; x <= _B_ + 0.0001; x += 0.1)
    {
        printf("%lf, ", ComputeSolution(x, _A_, _B_, n, uj));
    }

    printf("\n\n");

    printf("x = ");
    for (double x = _A_; x <= _B_ + 0.0001; x += 0.1)
    {
        printf("%lf, ", x);
    }

    free(uj);

    // printf("\n\nError on step h:\n");
    // const int MAX_NODES = 1000;
    // for (int m = 3; m < MAX_NODES; m++)
    // {
    //     double *u = new double[m + 1];
    //     GalerkinMethod(_A_, _B_, m, u);

    //     double *u_ex = new double[m];
    //     double *u_Gal = new double[m];
    //     for (int i = 0; i < m; i++)
    //     {
    //         //double x_cur = (ComputeRegularGridNode(_A_, _B_, i, m) + ComputeRegularGridNode(_A_, _B_, i + 1, m)) / 2;
    //         double x_cur = ComputeRegularGridNode(_A_, _B_, i, m);
    //         u_ex[i] = u_exact(x_cur);
    //         u_Gal[i] = ComputeSolution(x_cur, _A_, _B_, m, u);
    //     }
    //     double error_norm = ComputeErrorNorm(u_Gal, u_ex, m);
    //     //printf("error = %.20lf, h = %lf, m = %d\n", error_norm, (_B_ - _A_) / double(m), m);
    //     printf("%.15lf, ", error_norm);


    //     delete[] u_ex;
    //     delete[] u;
    //     delete[] u_Gal;
    // }

    // printf("\n\nStep:\n");

    // for (int m = 3; m < MAX_NODES; m++)
    //     printf("%.15lf, ", (_B_ - _A_) / (double)m);

    // return 0;
}

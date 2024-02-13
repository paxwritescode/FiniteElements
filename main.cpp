#include <iostream>
#include "GalerkinMethod.h"
#include "functions.h"
#include "test.h"

int main(void)
{
    for (double x = _A_; x <= _B_ + 0.0001; x += 0.1)
    {
        printf("%lf, ", GalerkinMethod(x, _A_, _B_, 10));
    }

    printf("\n\n");

    for (double x = _A_; x <= _B_ + 0.0001; x += 0.1)
    {
        printf("%lf, ", x);
    }

    printf("\n");

    Test();
    
    return 0;
}

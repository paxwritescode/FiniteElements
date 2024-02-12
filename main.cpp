#include <iostream>
#include "GalerkinMethod.h"
#include "functions.h"
#include "test.h"

int main(void)
{
    for (double x = _A_; x <= _B_; x += 0.1)
    {
        printf("%lf, ", GalerkinMethod(x, _A_, _B_, 5));
    }

    printf("\n");

    Test();
    
    return 0;
}

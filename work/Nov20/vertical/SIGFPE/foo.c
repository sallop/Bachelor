#include <stdio.h>

#define _GNU_SOURCE
#include <fenv.h>

int main(void)
{
        double x, y;

        int feenableexcept();
        feenableexcept(FE_ALL_EXCEPT);

        x= 0.0;
        y= 0.0;

        x= x/y;

        printf("%lf\n",x);
}

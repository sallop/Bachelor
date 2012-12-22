#include <stdio.h>

//#define _GNU_SOURCE
#include <fenv.h>
#include <signal.h>
#include <stdlib.h>

void fpehandler(int sig_num)
{
        signal(SIGFPE, fpehandler);
        printf("SIGFPE: floating point exception occured, exiting.\n");
        abort();
}

int main(void)
{
        double x, y;

     int feenableexcept();	// may be function prototype declare
 //  int printf(const char *str, ...);
        feenableexcept(FE_ALL_EXCEPT);
        signal(SIGFPE, fpehandler);

        x= 0.0;
        y= 0.0;

        x= x/y;

        printf("%lf\n",x);
}

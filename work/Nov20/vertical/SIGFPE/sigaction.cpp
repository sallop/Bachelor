#include <iostream>
#include <csignal>
#include <cmath>
#include <cerrno>
#include <fenv.h>

void fpe_handler(int signum, siginfo_t *info, void *ucontext)
{
  printf("%s\n", strsignal(signum));
  printf("signal emit at %p\n", info->si_addr);
  switch(info->si_code){
  case FPE_INTDIV:
    printf("\n");
    break;
  case FPE_INTOVF:
    printf("\n");
    break;
  case FPE_FLTDIV:
    printf("\n");
    break;
  case FPE_FLTOVF:
    printf("\n");
    break;
  case FPE_FLTUND:
    printf("\n");
    break;
  case FPE_FLTRES:
    printf("\n");
    break;
  case FPE_FLTINV:
    printf("\n");
    break;
  case FPE_FLTSUB:
    printf("\n");
    break;
  default:
    printf("?\n");
  }
  perror("hoge");
}

int main()
{
  double a;

  fexcept_t flag;
//  fegetexceptflag(&flag, FE_ALL_EXCEPT);
//  fegetexceptflag(&flag, FE_ALL_EXCEPT);

feenableexcept(FE_ALL_EXCEPT);
  //  fenv_t fenv;
  //  fegetenv(&fenv);
  //  fesetenv(&fenv);
  
  struct sigaction sa;
  //  sa.sa_flags = SA_NODEFER;  
  //  sa.sa_flags = SIGFPE;
  sa.sa_flags = SA_SIGINFO;
  sa.sa_handler = SIG_DFL;
  //sa.sa_mask = SIGFPE;
  sa.sa_sigaction = fpe_handler;
  //  sigaction(SIGFPE, &sa, NULL);
  //sigaction(SIGFPE, &sa, NULL);
  fetestexcept(FE_ALL_EXCEPT);
  //feraiseexcept(FE_DIVBYZERO|FE_INVALID);
  //raise(SIGSEGV);
  //raise(SIGFPE);

  //  feclearexcept(FE_ALL_EXCEPT);
  fetestexcept(FE_ALL_EXCEPT);  
  perror("before");

  sqrt(-1.2);
  fetestexcept(FE_ALL_EXCEPT);
  perror("after");  
  return 0;
}

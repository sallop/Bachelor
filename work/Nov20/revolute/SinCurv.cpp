#include <cstdio>
#include <cmath>
#include "SinCurv.hpp"
struct SinCurv;
		 
SinCurv::SinCurv(double _ttf,
		 double _th1_i,double _th1_f,
		 double _th1v_i, double _th1v_f)
  : ttf(_ttf),
    th1_i(_th1_i)  , th1_f(_th1_f),
    th1v_i(_th1v_i), th1v_f(_th1v_f)
{
  Amp1 = M_PI*(th1_f - th1_i)/(2.0*ttf) - M_PI*(th1v_f + th1v_i)/4.0;
  omg1 = M_PI/ttf;
}

double SinCurv::sin_th1(double t) const
{
  return th1_i + th1v_i*t
    + (th1v_f - th1v_i)*t*t/(2.0*ttf)
    + Amp1/omg1*(1.0 - cos(omg1*t));
}

double SinCurv::sin_th1v(double t) const
{
  return th1v_i + (th1v_f - th1v_i)*t/ttf + Amp1*sin(omg1*t);
}

double SinCurv::sin_th1a(double t) const
{
  return (th1v_f - th1v_i)/ttf + Amp1*omg1*cos(omg1*t);  
}

void SinCurv::print_value(FILE *_fp, double t)
{
  fprintf(_fp, "%lf %lf %lf %lf\n"
	  , t, sin_th1(t), sin_th1v(t), sin_th1a(t) );
}

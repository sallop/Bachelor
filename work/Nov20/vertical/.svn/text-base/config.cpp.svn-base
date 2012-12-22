#include "config.hpp"

//using namespace std;
Linear_curve_t::Linear_curve_t(double _a, double _b)
  : a(_a), b(_b)
{
  cout <<__FILE__<<"-"<<__FUNCTION__<<"-"<<__LINE__<< endl;
}

Sin_curve::Sin_curve(double _ttf,
		     double _th1_i,double _th1_f,
		     double _th1v_i, double _th1v_f)
  : ttf(_ttf),
    th1_i(_th1_i)  , th1_f(_th1_f),
    th1v_i(_th1v_i), th1v_f(_th1v_f)
{
  Amp1 = M_PI*(th1_f - th1_i)/(2.0*ttf) - M_PI*(th1v_f + th1v_i)/4.0;
  omg1 = M_PI/ttf;
}

double Sin_curve::sin_th1(double t) const
{
  return th1_i + th1v_i*t
    + (th1v_f - th1v_i)*t*t/(2.0*ttf)
    + Amp1/omg1*(1.0 - cos(omg1*t));
}

double Sin_curve::sin_th1v(double t) const
{
  return th1v_i + (th1v_f - th1v_i)*t/ttf + Amp1*sin(omg1*t);
}

double Sin_curve::sin_th1a(double t) const
{
  return (th1v_f - th1v_i)/ttf + Amp1*omg1*cos(omg1*t);  
}

void Sin_curve::print_value(FILE *_fp, double t)
{
  fprintf(_fp, "%lf %lf %lf %lf\n"
	  , t, sin_th1(t), sin_th1v(t), sin_th1a(t) );
}

BallScrew::BallScrew(double _J, double _m, double _pitch)
  : J(_J), m(_m), pitch(_pitch)
{
  frac = 2.0*M_PI/pitch;
  cout << "J=" << J << setw(2) 
       << "m=" << m << setw(2)
       << "f=" << frac << endl;
}

double BallScrew::getTau(double tha)
{
  return (J + m/(frac*frac))*tha;
}

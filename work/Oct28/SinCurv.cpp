#include <cstdio>
#include <cmath>
#include <iostream>
#include "SinCurv.hpp"

using namespace std;
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

SinCurv::SinCurv()
{
	ttf   = 0.0;
	th1_i = 0.0;
	th1_f = 0.0;
	th1v_i= 0.0;
	th1v_f= 0.0;
	Amp1  = 0.0;
	omg1  = 0.0;
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

SinCurv::SinCurv(const SinCurv &sc)
{
	std::cout << __FILE__ <<"-"<< __FUNCTION__ <<"-"<< __LINE__ << std::endl;
	std::cout << "No implementation this constructor.\n" << std::endl;
}

SinCurv& SinCurv::operator =(const SinCurv &sc)
{
	std::cout	<< __FILE__ <<"-"
				<< __FUNCTION__ <<"-"
				<< __LINE__ << std::endl;
	this->Amp1   = sc.Amp1  ;
	this->omg1   = sc.omg1  ;
	this->ttf    = sc.ttf   ;
	this->th1_i  = sc.th1_i ;
	this->th1_f  = sc.th1_f ;
	this->th1v_i = sc.th1v_i;
	this->th1v_f = sc.th1v_f;	
	return *this;
}

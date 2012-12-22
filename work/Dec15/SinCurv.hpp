#ifndef SINCURVE_H
#define SINCURVE_H
#include <iostream>
#include <cstdio>

using namespace std;
struct SinCurv{
  SinCurv(double _ttf,
	    double _th1_i , double _th1_f,
	    double _th1v_i, double _th1v_f);
  SinCurv();

  //  SinCurv(const SinCurv& sc);
  //  SinCurv& operator=(const SinCurv& sc); 

  double Amp1, omg1, ttf;
  double th1_i, th1_f, th1v_i, th1v_f;
  double sin_th1(double t) const;
  double sin_th1v(double t) const;
  double sin_th1a(double t) const;
  void print_value(FILE *_fp, double t);
};

#endif	// SINCURVE_H

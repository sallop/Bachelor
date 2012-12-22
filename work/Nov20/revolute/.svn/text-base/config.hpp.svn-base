#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <cmath>
#include <cstdio>
#include <iomanip>

using namespace std;

struct Linear_curve_t{
  Linear_curve_t(double _a, double _b);
  double a, b;
};

struct Conf_value{
  double x, y, z, vth, rth1, rth2;
};

struct Return_sim_value_t {double s, v, a;};

struct Sin_curve{
  Sin_curve(double _ttf,
	    double _th1_i , double _th1_f,
	    double _th1v_i, double _th1v_f);
  double Amp1, omg1, ttf;
  double th1_i, th1_f, th1v_i, th1v_f;
  double sin_th1(double t) const;
  double sin_th1v(double t) const;
  double sin_th1a(double t) const;
  void print_value(FILE *_fp, double t);
};

class BallScrew{
public:
  BallScrew(double _J, double _m, double _pitch);
  double J;// motor side of inerita
  double m;// mass of fuka side
  double frac;// fraction of gear teeth  = 2.0*M_PI/pitch;
  double pitch;// pitch of screw
  double getTau(double tha);// motor side of torque
};

#endif	// CONFIG_H

#ifndef __INCLUDE_VARIABLES_H__
#define __INCLUDE_VARIABLES_H__
#include <cmath>
#include "SinCurv.hpp"
#include "Manipulator.hpp"
struct Point{double x, y, z;}; // => mv hoge.h
struct Config { double ini, fin;}; // using th1[], th2[], th3[]};
struct Linear_curve_t { double a; double b; };
extern const double	g  ;
extern const double	dt ;
extern const double	Ra ;
extern const double	Im ;
extern const double	Dm ;
extern const double	Kt ;
extern const double	Kv ;
extern const double	HH ;
extern const double	RR ;
extern const double	LL ;
extern const double	L1 ;
extern const double	L2 ;
extern const double	L3 ;
extern const double	Lg1;
extern const double	Lg2;
extern const double	Lg3;
extern const double	IG1;
extern const double	IG2;
extern const double	m1 ;
extern const double	m2 ;
extern const double	m3 ;
extern const int	NN ;
extern const int	N3 ;
extern const double	Mf ;
extern const int    REPEAT ;

extern const double height;
extern const double radius;

extern const double b1;// b1=kv+Ra*Dm/Kt  0.052
extern const double b2;// b2=Ra*Im/Kt     0.000649
extern const double b3;// b3=Ra/Kt        76.08695

extern Point pnt_i ;
extern Point pnt_f ;
extern Point pnt_v0;
extern Point pnt_r0;
extern const double vth_i    ;
extern const double vth_f    ;
extern const double rth12    ;
extern const double line_i_r0;
extern const double line_i_v0;
extern const double rth1;
extern const Point pnt_r1;
//extern const  double rth2  = 2.0*acos(0.5*line_i_r0/L);
extern const double rth2   ;
extern const Config cfg_th1;
extern const Config cfg_th2;
extern const Config cfg_th3;
extern Linear_curve_t fk1, fk2;

extern double a1;
extern double a2;
extern double a3;
extern double a4;

extern ScaraRobot scara;

extern VRManipulator v1;// virtual manipulator
extern Revolute  r1;
extern Revolute  r2;// real manipulator 1, 2
extern Prismatic r3;// real manipulator 3
extern SinCurv  sc1;
extern SinCurv  sc2;
extern SinCurv  sc3;// initial trajectory
extern Energy e1;
extern Energy e2;
extern Energy e3;


extern double Ts;
extern double t1;

extern double ts[20];
extern double kiza1;
extern double kiza2;
extern double kiza3;
extern double Sth1[20];
extern double Sth2[20];
extern double Sth3[20];
extern double  th1[20];
extern double th1v[20];
extern double  th2[20];
extern double th2v[20];
extern double  th3[20];
extern double th3v[20];


#endif	// #define __INCLUDE_VARIABLES_H__

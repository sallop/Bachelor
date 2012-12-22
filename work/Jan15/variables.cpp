#include "variables.hpp"
#include "Manipulator.hpp"

const double	g   = 9.8;
//const double	dt  = 0.0001;// => 1000 <= (Ts / dt);
const double	dt  = 0.001;// => 1000 <= (Ts / dt);
//const double	dt  = 0.0001;// => 1000 <= (Ts / dt);

//const double	dt  = 0.01;// => 1000 <= (Ts / dt);
const double	Ra  = 3.5;
const double	Im  = 8.54e-06;
const double	Dm  = 7.89e-05;
const double	Kt  = 0.046;
const double	Kv  = 0.046;
const double	HH  = 0.08;
const double	RR  = 0.01;
const double	LL  = 0.08;
const double	L1  = 0.08;
const double	L2  = 0.08;
const double	L3  = 0.08;
const double	Lg1 = 0.04;
const double	Lg2 = 0.04;
const double	Lg3 = 0.04;
const double	IG1 = 1.73e-05;
const double	IG2 = 1.73e-05;
const double	m1  = 0.0202;
const double	m2  = 0.0202;
const double	m3  = 0.0202;	//#define m3	0.202//#define m3	2.020
const int	NN  = 7;	//#define NN	7
const int       N3  = static_cast<int>(NN * 0.5 + 0.5); // center value
const double	Mf  = 0.013;	// 0.03 coulomb friction
//const int    REPEAT = 60;	// 80, 101(limit), 102(illeagle),
const int    REPEAT = 10;	// 80, 101(limit), 102(illeagle), 

const double b1 = Kv+Ra*Dm/Kt;// b1=Kv+Ra*Dm/Kt  0.052
const double b2 = Ra*Im/Kt   ;// b2=Ra*Im/Kt     0.000649
const double b3 = Ra/Kt      ;// b3=Ra/Kt        76.08695

//Manipulator init position configuration
//const double height = 0.0618;
const double height = 0.05;

const double radius = 0.04  ;

// Point pnt_i  = { 0.05        , 0.10       , 0.0  };// real end_effector start
// Point pnt_f  = {-0.05        , 0.10       , 0.0  };// real end_effector end
// Point pnt_v0 = {-0.05 - 0.005, 0.10 + 0.05, 0.0  };// virtual base point
// Point pnt_r0 = { 0.00        , 0.00       , height};// real base point
const Point pnt_i(  0.05        , 0.10       , 0.0);// real end_effector start
const Point pnt_f( -0.05        , 0.10       , 0.0);// real end_effector end
const Point pnt_c(  0.00        , 0.05       , 0.0); // E1, E2, E3 change point.
const Point pnt_v0(-0.05 - 0.005, 0.10 + 0.05, 0.0);// virtual base point
//const Point pnt_v0(0.0        , 0.10       , 0.0);// virtual base point
const Point pnt_r0( 0.00        , 0.00       , height);// real base point

const double vth_i = atan2(pnt_i.y - pnt_v0.y, pnt_i.x - pnt_v0.x);
const double vth_f = atan2(pnt_f.y - pnt_v0.y, pnt_f.x - pnt_v0.x);
const double rth12 = atan2(pnt_i.y - pnt_r0.y, pnt_i.x - pnt_r0.x);
const double line_i_r0 = hypot(pnt_i.x - pnt_r0.x, pnt_i.y - pnt_r0.y);
const double line_i_v0 = hypot(pnt_i.x - pnt_v0.x, pnt_i.y - pnt_v0.y);  

const double rth1  = rth12 - acos(0.5*line_i_r0/L1);
// const Point pnt_r1 = {
//   pnt_r0.x + L1*cos(rth1), pnt_r0.y + L1*sin(rth1), pnt_r0.z,    
//};
const Point pnt_r1( pnt_r0.x + L1*cos(rth1),
		    pnt_r0.y + L1*sin(rth1),
		    pnt_r0.z );

//const  double rth2  = 2.0*acos(0.5*line_i_r0/L);
const double rth2 = atan2(pnt_i.y - pnt_r1.y, pnt_i.x - pnt_r1.x);

const struct Config cfg_th1 = { 0.0   , height }; // distance
const struct Config cfg_th2 = { vth_i , vth_f  }; // degree
const struct Config cfg_th3 = { height, 0.0    }; // distance

// variable
struct Linear_curve_t fk1, fk2, fk3, fk4;
double a1;
double a2;
double a3;
double a4;

ScaraRobot scara;

VRManipulator v1;// virtual manipulator
Revolute  r1, r2;// real manipulator 1, 2
Prismatic r3    ;// real manipulator 3
Energy e1;	 // r1
Energy e2;	 // r2
Energy e3;	 // r3

// initial trajectory
SinCurv sc1; // Up
SinCurv sc2; // Direct
SinCurv sc3; // Down

double Ts;// Working time
double t1;// unit

double ts[20];// Time separatoar
double kiza1; // IDP kizami haba. r1.
double kiza2; // r2
double kiza3; // r3
double Sth1[20];		// Saiteki kidou r1.th
double Sth2[20];		// r2.th
double Sth3[20];		// r3.z0
double  th1[20];
double th1v[20];
double  th2[20];
double th2v[20];
double  th3[20];
double th3v[20];

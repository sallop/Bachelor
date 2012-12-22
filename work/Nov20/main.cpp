#include <cstdio>
#include <cfloat>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>
#include <fenv.h>
#include "GraphicalPrimitive.hpp"
#include "Manipulator.hpp"
#include "SinCurv.hpp"
#include "Energy.hpp"
//#include "ReadTrajectory.hpp"
//#include "debug.h"

using namespace std;

#define g	9.8
//#define g	0.0
#define dt	0.001
#define Ra	3.5
#define Im	8.54e-06
#define Dm	7.89e-05
#define Kt	0.046
#define Kv	0.046
#define HH	0.08
#define RR	0.01
#define LL	0.08
#define L1	0.08
#define L2	0.08
#define L3	0.08
#define Lg1	0.04
#define Lg2	0.04
#define Lg3	0.04
#define IG1	1.73e-05
#define IG2	1.73e-05
#define m1	0.0202
#define m2	0.0202
#define m3	0.0202
//#define m3	0.202
//#define m3	2.020
#define NN	7    // center value

//#define Mf      0.03 // not viscosity resistance. This is coulomb friction
const double Mf = 0.013;
const int REPEAT = 60;
//const int REPEAT = 80;
//const int REPEAT = 101;	// limit
//const int REPEAT = 102;	// ill
int plot_flags[1+REPEAT] = {1};//Animation flag. offset is 1

// configuration variable
struct Linear_curve_t { double a; double b; } fk1, fk2;
struct Variable_t { double s; double v; double a; };
// s is initialize of `self`.

double th1_i, th1_f, th1v_i, th1v_f;
double th2_i, th2_f, th2v_i, th2v_f;
double th3_i, th3_f, th3v_i, th3v_f;

// global values
double t;			// time
Energy e1, e2, e3;// using manipulator r1, r2, r3{ev, ia, tau, energy};
// number is joint postion.
double Ene123;
double Enemin1[20], Enemin2[20], Enemin3[20];
VRManipulator v1;		// virtual manipulator
Revolute r1, r2;		// real manipulator 1, 2
Prismatic r3;			// real manipulator 3
SinCurv sc1, sc2, sc3;		// initial trajectory
// using DP_variable
int Mat[20][20][20], N3;
 // armature cofficient
double b1;// b1=kv+Ra*Dm/Kt
double b2;// b2=Ra*Im/Kt
double b3;// b3=Ra/Kt
double Energ[10][10][10], enemin1, enemin2, Enemin123[20];
double Sth1[20], th1[20], th1v[20], th1a[20];
double Sth2[20], th2[20], th2v[20], th2a[20];
double Sth3[20], th3[20], th3v[20], th3a[20];
double ts[20], Ts, t1;
int ntmax[1+4]; // static_cast<int>(ts[]/dt);
// ts[offset], offset=1, Ts=whole time, t1=one part of times
double a1, a2, a3, a4;
double kiza1, kiza2, kiza3;

int sat[1 + 6], at[1 + 6];
int nnn; // repeat times

// using DP_function
void repeat ();
void tansaku();
void optimal_result(int make_animation_flag);
void make_gif_file();
void settei (int step);
void sim(int step);
void sim2(int step);
void sim3(int step);

double tha_1, thv_1, th_1;	// return sim() parametar
double tha_2, thv_2, th_2;	// return sim() parametar
double tha_3, thv_3, th_3;	// return sim() parametar
void Hikaku1 (int j);
void zu(ofstream &ofs);
void sinc();

inline void
graph_volt_constituent(ofstream& ofs, double t, double ev,
		       double term1, double term2, double term3)
{
  ofs << t     <<" "
      << ev    <<" "
      << term1 <<" "
      << term2 <<" "
      << term3 <<" "
      << endl;
}

inline void
graph_standard(ofstream& ofs, double t,
	       double z, double zv, double za,
	       double ev, double ia, double tau,
	       double energy, double ev_ia)
{
  ofs << t      << " "		// 1
      << z      << " "		// 2
      << zv     << " "		// 3
      << za     << " "		// 4
      << ev     << " "		// 5
      << ia     << " "		// 6
      << tau    << " "		// 7
      << energy << " "		// 8
      << ev_ia  << " "		// 9
      << std::endl;
}
void
reinitialize()
{
  t = 0.0;
 
  e1 = Energy();
  r1.th = 0.0; r1.thv = 0.0; r1.tha = 0.0;

  e2 = Energy();
  r2.th = 0.0; r2.thv = 0.0; r2.tha = 0.0;
   
  e3 = Energy();
  r3.z0 = 0.0; r3.zv = 0.0; r3.za = 0.0;
}

void
init()
{
  // \< Don't change follow lines sequence.
  //  const double height = 0.0618;

  // after modify {const double height = 0.0618*2; => 0.0618*5}
  const double height = 0.0618;
  //  const double height = 0.0618*2;
  //  const double height = 0.0618*5;
  //  const double height = 0.0618*5;

  // optimal pulley's radius = `r = 0.04` when rouund-trip.
  const double radius = 0.04;

  const struct Point{
    double x, y, z;
  } pnt_i  = { 0.05        , 0.10       , 0.0   } // real end_effector start
  , pnt_f  = {-0.05        , 0.10       , 0.0   } // real end_effector end
  , pnt_v0 = {-0.05 - 0.005, 0.10 + 0.05, 0.0   } // virtual base point
  , pnt_r0 = { 0.00        , 0.00       , height}; // real base point
  
  const double
    vth_i     = atan2(pnt_i.y - pnt_v0.y, pnt_i.x - pnt_v0.x),
    vth_f     = atan2(pnt_f.y - pnt_v0.y, pnt_f.x - pnt_v0.x),
    rth12     = atan2(pnt_i.y - pnt_r0.y, pnt_i.x - pnt_r0.x),
    line_i_r0 = hypot(pnt_i.x - pnt_r0.x, pnt_i.y - pnt_r0.y),
    line_i_v0 = hypot(pnt_i.x - pnt_v0.x, pnt_i.y - pnt_v0.y);  

  const double rth1  = rth12 - acos(0.5*line_i_r0/L1);
  const struct Point pnt_r1 = {
    pnt_r0.x + L1*cos(rth1), pnt_r0.y + L1*sin(rth1), pnt_r0.z,    
  };
  
  //  double rth2  = 2.0*acos(0.5*line_i_r0/L);
  const double rth2 = atan2(pnt_i.y - pnt_r1.y, pnt_i.x - pnt_r1.x);

  const struct Config {
    double ini, fin;// using th1[], th2[], th3[]
  } cnf_th1 = { 0.0   , height } // distance
  , cnf_th2 = { vth_i , vth_f  } // degree
  , cnf_th3 = { height, 0.0    }; // distance
  // -- Don't change sequence \>
    
  // initialize global value
  //  t  = 0.0; Ts = 0.64; t1 = Ts/8.0;// pull pulley_test2.cpp
  //  t  = 0.0; Ts = 0.40; t1 = Ts/8.0;// pull pulley_test2.cpp
  //  t  = 0.0; Ts = 0.16; t1 = Ts/8.0;// pull pulley_test2.cpp
  t  = 0.0; Ts = 0.12; t1 = Ts/8.0;// pull pulley_test2.cpp
  //  t  = 0.0; Ts = 0.64*2.0; t1 = Ts/8.0;// pull pulley_test2.cpp
  //  t  = 0.0; Ts = 0.1; t1 = Ts/8.0;// pull pulley_test2.cpp
  //  t  = 0.0; Ts = 0.16; t1 = Ts/8.0;// pull pulley_test2.cpp  
  //  t1 = Ts/(double)(DIV)
  //  for(int i=0; i <= DIV; i++){ ts[i] = t1*i; }
  for(int i=0; i <= 8; i++){ ts[i] = t1*i; }
  
  b1 = Kv + Ra*Dm/Kt;
  b2 = Ra*Im/Kt;
  b3 = Ra/Kt;
  
  //  N3 = static_cast<int>(NN * 0.5 + 0.5);
  N3 = static_cast<int>(NN * 0.5 + 0.5);

  //height= 0.0618

  // kiza1 = 0.005;
  // kiza1 = 0.0001;
  // kiza1 = 0.0002;
  // kiza1 = 0.0003;
  // kiza1 = 0.0004;
  // kiza1 = 0.1;
  kiza1 = 0.01;
  
  //  kiza1   = 0.001;
  kiza2 = 0.01;
  kiza3 = 0.01;
  
  // Revolute(mass, length, x0, y0, z0, rth);
  v1 = VRManipulator( 0.0, line_i_v0,
		      pnt_v0.x, pnt_v0.y, pnt_v0.z,
		      vth_i);
  //  v1.setXYZ(pnt_i.x, pnt_i.y, pnt_i.z);
  v1.setXYZ(pnt_i.x, pnt_i.y, pnt_i.z);  

  // pnt_r0 is a joint base position
  r1 = Revolute(m1, L1, pnt_r0.x, pnt_r0.y, pnt_r0.z, rth1);
  r1.setXYZ(pnt_r1.x, pnt_r1.y, pnt_r1.z);
  
  r2 = Revolute(m2, L2, pnt_r1.x, pnt_r1.y, pnt_r1.z, rth2);
  r2.setXYZ(pnt_i.x, pnt_i.y, pnt_i.z);
  
  r3 = Prismatic(m3, L3, pnt_i.x, pnt_i.y, pnt_i.z, radius);
  r3.setXYZ(pnt_i.x, pnt_i.y, pnt_i.z - r3.l);
  // x0, y0, z0);	// assignment: 仮想マニピュレータ初期化
  // 代入。 not use copy constructor
  //  sc1 = SinCurv(0.64, 0.0, 1.0, 0.0, 0.0);// 上昇 pull lab/work2.cpp
  //  sc1 = SinCurv(0.64, 0.0, 0.0618, 0.0, 0.0);// 上昇 pull lab/work2.cpp
  //  sc3 = SinCurv(0.64, 1.0,    0.0, 0.0, 0.0);// 下降
  
  sc1 = SinCurv(  Ts, cnf_th1.ini, cnf_th1.fin, 0.0, 0.0);// up
  sc2 = SinCurv(  Ts, cnf_th2.ini, cnf_th2.fin, 0.0, 0.0); // prismatic
  sc3 = SinCurv(  Ts, cnf_th3.ini, cnf_th3.fin, 0.0, 0.0);// down

  th_1  = cnf_th1.ini;		// initialize value
  thv_1 = sc1.sin_th1v(0.0);	// 
  tha_1 = sc1.sin_th1a(0.0);	// 

  th_2  = cnf_th2.ini;		// initialize value
  thv_2 = sc2.sin_th1v(0.0);	// 
  tha_2 = sc2.sin_th1a(0.0);	// 

  th_3  = cnf_th3.ini;		// initialize value
  thv_3 = sc3.sin_th1v(0.0);	// 
  tha_3 = sc3.sin_th1a(0.0);	// 
  {// Setting Sth1[]
    ofstream ofs("dir-dat/sin_th1.dat");
    for(int i=0; i <= 8; ++i){//Magick
      Sth1[i] = sc1.sin_th1(ts[i]);
      Sth2[i] = sc2.sin_th1(ts[i]);
      Sth3[i] = sc3.sin_th1(ts[i]);
      ofs << ts[i] <<" "<< Sth1[i] << std::endl;
    }
    ofs.close();
  }

  // Nov 18
  // ReadTrajectory rt("trajectory.dat");
  // set value in Sth[i] (or Rth[i]) the interpolated data
  
  ntmax[0] = 0;
  ntmax[1] = static_cast<int>(ts[3]/dt);
  ntmax[2] = static_cast<int>(ts[4]/dt);
  ntmax[3] = static_cast<int>(ts[5]/dt);
  ntmax[4] = static_cast<int>(ts[8]/dt);//Magick
  
  {// initialize theta
    //    th1_i  = vth_i; th1_f  = vth_f ;
    th1_i = cnf_th1.ini; th1_f = cnf_th1.fin;
    th1v_i = 0.0       ; th1v_f = 0.0       ;
    th1[0] = th1_i     ; th1[8] = th1_f;// th1[DIV] = th1_f;
    th1v[0] = 0.0      ; th1v[8] = th1v_f; //    th1v[DIV] = th1v_f;
    // th2[]
    th2_i = cnf_th2.ini; th2_f = cnf_th2.fin;
    th2v_i = 0.0       ; th2v_f = 0.0   ;
    th2[0] = th2_i     ; th2[8] = th2_f;// th2[DIV] = th2_f;
    th2v[0] = 0.0      ; th2v[8] = th2v_f; //    th2v[DIV] = th2v_f;
    // th3[]
    th3_i = cnf_th3.ini; th3_f = cnf_th3.fin;
    th3v_i = 0.0       ; th3v_f = 0.0   ;
    th3[0] = th3_i     ; th3[8] = th3_f; // th3[DIV] = th3_f;
    th3v[0] = 0.0      ; th3v[8] = th3v_f;  // th3v[DIV] = th3v_f;
  }

  {//init coefficient
    a1 = r1.Ig + r1.m*r1.r*r1.r + r2.m*r1.l*r1.l + r3.m*r1.l*r1.l;
    a2 = r2.Ig + r2.m*r2.r*r2.r + r3.m*r2.l*r2.l;
    a3 = r3.Iz;
    a4 = r2.m*r1.l*r2.r + r3.m*r1.l*r2.l;   
  }
  
  ofstream ofs("init_position.dat");
  zu(ofs);
}

void
view_animation( std::string prefix )
{
  FILE *gp = popen("gnuplot","w");
  fprintf(gp, "set xrange[%lf:%lf]\n", -0.08, 0.10);
  fprintf(gp, "set yrange[%lf:%lf]\n", -0.02, 0.16);
  fprintf(gp, "set zrange[%lf:%lf]\n",  0.00, 0.20);
  //  fprintf(gp, "set view %lf, %lf\n", 67.0, 133.0);
  fprintf(gp, "set view %lf, %lf\n", 0.0, 0.0);
  fprintf(gp, "set terminal gif\n");
  
  
  for(int ii=10; ii < static_cast<int>(ts[8]/dt); ii+=10)
    {//Magick
      fprintf(gp, "set output 'dir-gif/%s-%05d.gif\n"
	      , prefix.c_str(), ii);
    
      fprintf(gp, "splot 'dir-dat/%s-%05d.dat' w l\n"
	      , prefix.c_str(), ii);
      
      fprintf(gp, "pause 0.01\n");      
    }

  pclose(gp);
}

void
zu()
{
  cout <<__FILE__<<"_"<<__FUNCTION__<<"_"<<__LINE__<< endl;
}

void
zu(std::ofstream& ofs)
{
  cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<< endl;
  Circle3D real1_p0(r1.x0, r1.y0, r1.z0, 0.01);
  Circle3D real1_p1(r1.x , r1.y , r1.z0, 0.01);

  Circle3D real2_p0(r2.x0, r2.y0, r2.z0, 0.01);
  Circle3D real2_p1(r2.x , r2.y , r2.z0, 0.01);
  Circle3D real2_p2(r3.x0, r3.y0, r3.z0, 0.01);

  Circle3D real3_p0(r3.x0, r3.y0, r3.z0, 0.01);
  Circle3D real3_p1(r3.x0, r3.y0, r3.z , 0.01);
  
  Line3D real1_ln(r1.x0, r1.y0, r1.z0, r1.x,  r1.y , r1.z0 );
  Line3D real2_ln(r2.x0, r2.y0, r1.z0, r2.x,  r2.y , r1.z0 );
  Line3D real3_ln(r3.x0, r3.y0, r3.z0, r3.x0, r3.y0, (r3.z = r3.z0 + r3.l));
  
  Circle3D virt1_p0(v1.x0, v1.y0, v1.z0, 0.01);
  Circle3D virt1_p1(v1.x , v1.y , v1.z0, 0.01);

  Line3D virt1_ln(v1.x0, v1.y0, v1.z0,
		  v1.x , v1.y , v1.z0);
  
  // r1
  real1_p0.draw(ofs);
  real1_p1.draw(ofs);// real1_p2.draw(ofs);
  real1_ln.draw(ofs);
  // r2
  real2_p0.draw(ofs);
  real2_p1.draw(ofs);// real2_p2.draw(ofs);
  real2_ln.draw(ofs);
  // r3
  real3_p0.draw(ofs);
  real3_p1.draw(ofs);// real3_p2.draw(ofs);
  real3_ln.draw(ofs);
  // v1
  virt1_p0.draw(ofs);
  virt1_p1.draw(ofs);// virt1_p2.draw(ofs);
  virt1_ln.draw(ofs);  
}

void
En_up(Prismatic& r3)
{
  e3.tau = r3.m*(r3.za + g)*r3.radius;
  // pulley's radius and torque.
  e3.ev  = b1*(r3.zv/r3.radius) + b2*(r3.za/r3.radius) + b3*e3.tau;
  // Nov18 19:00
  // follow paragraph is comment out 
  // Nov19
  //if (r3.zv > 0.0) {
  //  e3.ev = e3.ev + b3*Mf;
  //} else if (r3.zv < 0.0) {
  //  e3.ev = e3.ev - b3*Mf;
  //}
  e3.ia  = (e3.ev - Kv*r3.zv/r3.radius)/Ra;
  
  if (e3.ev*e3.ia > 0.0)
    {
      e3.energy = e3.energy + e3.ev*e3.ia*dt;
    }
}

void
En_down(Prismatic& r3)
{
  double F = r3.m*(r3.za + g);
  e3.tau   = F*r3.radius;  	// pulley's radius and torque.
  e3.ev    = b1*(r3.zv/r3.radius) + b2*(r3.za/r3.radius) + b3*e3.tau;
  e3.ia    = (e3.ev - Kv * r3.zv/r3.radius)/Ra;

  //if (r3.zv > 0.0){
  //    e3.ev = e3.ev + b3*Mf;
  //}
  //else if (r3.zv < 0.0){
  //  e3.ev = e3.ev - b3*Mf;
  //}

  if (e3.ev*e3.ia > 0.0)
    {
      e3.energy = e3.energy + e3.ev*e3.ia*dt;
      //Nov5 11:20    printf("line=%d, e3.energy=%lf\n", __LINE__, e3.energy);
    }
}

double 
En_prism(const struct Return_sim_value_t &th_)
{
}

/*
 * Nov 21
 * modify En_prism(const struct Return_sim_value_t &th_){}
 */
void
En_direct(Revolute& r1, Revolute& r2, VRManipulator& v1)
{
  std::cout << "Nothing do this function. "
	    << __FUNCTION__
	    << std::endl;
  register double detJ, J11, J12, J21, J22, coef;
  register double cth1, cth12, cth2;
  register double sth1, sth12, sth2;  

  cth1 = cos(r1.th); cth12 = cos(r1.th + r2.th); cth2 = cos(r2.th);
  sth1 = sin(r1.th); sth12 = sin(r1.th + r2.th); sth2 = sin(r2.th);

  fk1.a = (fin.y - ini.y)/(fin.x - ini.x);
  fk1.b = ini.y;
  fk2.a = tan(th_.s);
  fk2.b = v1.y0 - fk2.a*v1.x0 ;
  fxx = (fk2.b - fk1.b)/(fk1.a - fk2.a);
  fyy = fk2.a*fxx + fk2.b;

  v1.l = hypot(fxx-v1.x0 , fyy-v1.y0);
  coef = (fk1.a*tan(th_.s) + 1.0)/(fk1.a - tan(th_.s));
  v1.lv = coef*v1.l*th_.v;
  v1.xv1 = v1.lv*cos(th_.s) - v1.l*th_.v*sin(th_.s);
  v1.yv1 = v1.lv*sin(th_.s) + v1.l*th_.v*cos(th_.s);

  keirosokudo = sqrt(v1.xv1*v1.xv1 + v1.yv1*v1.yv1);
  J11 = -r1.l*sth1 - r2.l*sth12;
  J12 = -r2.l*sth12;
  J21 =  r1.l*cth1 + r2.l*cth12;
  J22 =  r2.l*cth12;

  detJ = det(J11,J12,J21,J22);
  r1.thv = det(v1.xv1, J12, v1.yv1, J22)/detJ;
  r2.thv = det(   J11, v1.xv1, J21, v1.yv1)/detJ;

  PL3 = sqrt(fxx*fxx + fyy*fyy);

  r1.th = atan2(fyy,fxx) - acos(0.5*PL3/r1.l);
  r2.th = 2.0*acos(0.5*PL3/r1.l);

  v1.la = coef*(v1.l*th_.a + 2.0*v1.lv*th_.v) + v1.l*th_.v*th_.v;

  v1.xa1 = v1.la*cos(th_.s) - 2.0*v1.lv*th_.v*sin(th_.s)
    - v1.l*th_.a*sin(th_.s) - v1.l*th_.v*th_.v*cos(th_.s);
  v1.ya1 = v1.la*sin(th_.s) + 2.0*v1.lv*th_.v*cos(th_.s)
    + v1.l*th_.a*cos(th_.s) - v1.l*th_.v*th_.v*sin(th_.s);

  XXa = v1.xa1 + r1.thv*r1.thv*r1.l*cth1
    + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*cth12;
  YYa = v1.ya1 + r1.thv*r1.thv*r1.l*sth1
    + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*sth12;

  r1.tha = det(XXa, J12, YYa, J22)/detJ;
  r2.tha = det(J11, XXa, J21, YYa)/detJ;

  cth1 = cos(r1.th); cth12 = cos(r1.th + r2.th); cth2 = cos(r2.th);
  sth1 = sin(r1.th); sth12 = sin(r1.th + r2.th); sth2 = sin(r2.th);

  //tau_1 = (a1+a2+a3+2.0*a4*cth2)*r1.tha + (a2+a3+a4*cth2)*r2.tha
  //  - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*sth2;
  //tau_2 = (a2 + a3 + a4*cth2)*r1.tha
  //  + (a2 + a3)*r2.tha + a4*r1.thv*r1.thv*sth2;
	e1.tau = (a1+a2+a3+2.0*a4*cth2)*r1.tha + (a2+a3+a4*cth2)*r2.tha 
	- a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*sth2;
	e2.tau =(a2 + a3 + a4*cth2)*r1.tha + (a2 + a3)*r2.tha + a4*r1.thv*r1.thv*sth2;
	
  //ev1 = b1*r1.thv + b2*r1.tha + b3*tau_1; ia1 = (ev1 - Kv*r1.thv)/Ra;
  //ev2 = b1*r2.thv + b2*r2.tha + b3*tau_2; ia2 = (ev2 - Kv*r2.thv)/Ra;
e1.ev = b1*r1.thv + b2*r1.tha + b3*tau_1; e1.ia = (e1.ev - Kv*r1.thv)/Ra;
e2.ev = b1*r2.thv + b2*r2.tha + b3*tau_2; e2.ia = (e2.ev - Kv*r2.thv)/Ra;
  //if (ev1*ia1 > 0.0){ Ene_1 = Ene_1 + ev1*ia1*dt; }
  //if (ev2*ia2 > 0.0){ Ene_2 = Ene_2 + ev2*ia2*dt; }

  //Ene123 = Ene_1 + Ene_2 + Ene_3;
  //Ene123 = Ene_1 + Ene_2;
Ene122 = e1.energy + e2.energy + e3.energy;
  return Ene123;//Ene123;  
}

void
En_sim(int ntmax, int cur)
{
  reinitialize();
  for(int nt = 1; nt <= ntmax; nt++)
    {
      t += dt;
      sim(cur);// set th_1, thv_1, tha_1
      r3.z0 = th_1 ;
      r3.zv = thv_1;
      r3.za = tha_1; 

      // if (r3.zv < 0.0)
      // {
      // 	  e3.energy = DBL_MAX;
      // 	  break;
      // }

      En_up(r3);
      // if (nt < )     En_direct(r1, r2, v1);
      // if (nt < )     En_down(r3);
    }
  //  Energ[ at[cur] ][ at[cur+1] ][ at[cur+2] ]
  // = Ene123 = e3.energy;
  Energ[ at[cur] ][ at[cur+1] ][ at[cur+2] ] = e3.energy;
}


int
main(int argc, char *argv[])
{
  //  feenableexcept (FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  
  
  Manipulator m(m1, LL, 0.0, 0.0, 0.0);
  // Manipulator(mass, length, x0, y0, z0);
  init();
  sinc();
  view_animation("e1");		// view sinc animation


  assert(false);
  
  {// plot enemin2 block
    ofstream ofs("dir-dat/enemin2.dat");
    //    int plot_flags[1+60] = {0};
    //Nov18    int plot_flags[1+60] = {1};
    int plot_flags[1+REPEAT] = {1};
    for(int ii=0; ii < REPEAT; ++ii)
      plot_flags[ii] = 1;
    
    plot_flags[ 1] = 1; plot_flags[10] = 1; plot_flags[20] = 1;
    plot_flags[30] = 1; plot_flags[39] = 1; plot_flags[40] = 1;
    plot_flags[50] = 1; plot_flags[60] = 1;
    //for (int nnn = 1; nnn <= 13; nnn++) {
    //    for (nnn = 1; nnn <= 60; nnn++)
    for (nnn = 1; nnn <= REPEAT; nnn++)
      {//Magick
	printf("make animation flag=%d\n", plot_flags[nnn]);
	repeat();
	enemin2 = 0.0;
	{// tansaku minimum energy
	  tansaku();
	  optimal_result( plot_flags[nnn] );
	  //if( plot_flags[nnn] )
	  //  {
	  //    optimal_result(1);
	  //  }
	  //else
	  //  {
	  //    optimal_result(0);
	  //  }
	}
      
	Enemin2[nnn] = enemin2;

	if ( (nnn > 15) && (Enemin2[nnn] > Enemin2[nnn - 1]) )
	  {//Magick
	    kiza1 *= 1.2; //Magick
	  }

	ofs << nnn <<" "<< enemin2 <<" "<< kiza1 << endl;
      }
  }//end plot enemin2

  repeat ();
  enemin2 = 0.0;
  {
    tansaku();
    optimal_result(1);
  }
  // if speed is more precious replace this function `make_animation`.
  make_gif_file();

  
  view_animation("dp");
  
  return 0;
}

void
sinc()//Magick
{
  int fnum;
  double time_stack[4] = {0.0, Ts*1, Ts*2, Ts*3 };
  ofstream ofs, ofs2;
  ofs.open("dir-dat/sinc.dat");
  ofs2.open("dir-dat/sinc_ev.dat");

  reinitialize();
  for (int i = 0, fno = 0; i < static_cast <int>(ts[8]/dt); ++i, ++fnum)
    {//Magic
      t += dt;
      th_1  = sc1.sin_th1(t);
      thv_1 = sc1.sin_th1v(t);
      tha_1 = sc1.sin_th1a(t);

      r3.z0 = th_1 ;
      r3.zv = thv_1;
      r3.za = tha_1;
   
      En_up(r3);

      if (i % 10 == 0)
	{
	  char fname[32];
	  //	  sprintf(fname, "dir-dat/e1-%05d.dat", i);
	  sprintf(fname, "dir-dat/e1-%05d.dat", fno);
	  ofstream fout(fname);
	  zu(fout);
	}

      graph_standard(ofs, t,
		     th_1, thv_1, tha_1,
		     e3.ev, e3.ia, e3.tau,
		     e3.energy, e3.ev*e3.ia);
      
      graph_volt_constituent(ofs2, t, e3.ev,
			     b1*(r3.zv/r3.radius),
			     b2*(r3.za/r3.radius),
			     b3*(e3.tau));
    }

  for (int i = 0; i < static_cast <int>(ts[8]/dt); ++i, ++fnum)
    {//Magic
      t += dt;
      th_1  = sc1.sin_th1(t);
      thv_1 = sc1.sin_th1v(t);
      tha_1 = sc1.sin_th1a(t);

      v1.th  = th_1 ;
      v1.thv = thv_1;
      v1.tha = tha_1;
      
      En_direct(r1, r2, v1);

      if (i % 10 == 0)
	{
	  char fname[32];
	  //  sprintf(fname, "dir-dat/e1-%05d.dat", i);
	  sprintf(fname, "dir-dat/e1-%05d.dat", fnum);
	  ofstream fout(fname);
	  zu(fout);
	}

      graph_standard(ofs, t,
		     th_1, thv_1, tha_1,
		     e3.ev, e3.ia, e3.tau,
		     e3.energy, e3.ev*e3.ia);
      
      graph_volt_constituent(ofs2, t, e3.ev,
			     b1*(r3.zv/r3.radius),
			     b2*(r3.za/r3.radius),
			     b3*(e3.tau));
    }

  for (int i = 0; i < static_cast <int>(ts[8]/dt); ++i, ++fnum)
    {//Magic
      t += dt;
      th_1  = sc1.sin_th1(t);
      thv_1 = sc1.sin_th1v(t);
      tha_1 = sc1.sin_th1a(t);

      r3.z0 = th_1 ;
      r3.zv = thv_1;
      r3.za = tha_1;
   
      En_down(r3);

      if (i % 10 == 0)
	{
	  char fname[32];
	  sprintf(fname, "dir-dat/e1-%05d.dat", i);
	  ofstream fout(fname);
	  zu(fout);
	}

      graph_standard(ofs, t,
		     th_1, thv_1, tha_1,
		     e3.ev, e3.ia, e3.tau,
		     e3.energy, e3.ev*e3.ia);
      
      graph_volt_constituent(ofs2, t, e3.ev,
			     b1*(r3.zv/r3.radius),
			     b2*(r3.za/r3.radius),
			     b3*(e3.tau));
    }

  
  ofs.close();
  ofs2.close();
}

inline void
sim(int step)//Magick
{
  if (t < ts[1])
    {
      tha_1 = th1a[1];
      thv_1 = tha_1 *	t;
      th_1  = 0.5 * tha_1 * t * t + th1_i;
    }
  if (ts[1] < t && t < ts[2])
    {
      tha_1 = th1a[2];
      thv_1 = th1v[1] + tha_1 * (t - ts[1]);
      th_1  = th1[1] + th1v[1]*(t - ts[1]) + 0.5*tha_1*pow((t - ts[1]), 2);
    }
  if (ts[2] < t && t < ts[3])
    {
      tha_1 = th1a[3];
      thv_1 = th1v[2] + tha_1 * (t - ts[2]);
      th_1  = th1[2] + th1v[2]*(t - ts[2]) + 0.5*tha_1*pow((t - ts[2]), 2);
    }

  if (1 == step)
    return;

  if (ts[3] < t && t < ts[4])
    {
      tha_1 = th1a[4];
      thv_1 = th1v[3] + tha_1*(t - ts[3]);
      th_1  = th1[3] + th1v[3]*(t - ts[3]) + 0.5 * tha_1 * pow ((t - ts[3]), 2);
    }
  if (2 == step)
    return;

  if (ts[4] < t && t < ts[5])
    {
      tha_1 = th1a[5];
      thv_1 = th1v[4] + tha_1 * (t - ts[4]);
      th_1  = th1[4] + th1v[4] * (t - ts[4])
	+ 0.5 * tha_1 *	pow ((t - ts[4]), 2);
    }
  if (3 == step)
    return;

  if (ts[5] < t && t < ts[6])
    {
      tha_1 = th1a[6];
      thv_1 = th1v[5] + tha_1 * (t - ts[5]);
      th_1  = th1[5] + th1v[5]*(t - ts[5]) + 0.5*tha_1*pow((t - ts[5]), 2);
    }
  if (ts[6] < t && t < ts[7])
    {
      tha_1 = th1a[7];
      thv_1 = th1v[6] + tha_1 * (t - ts[6]);
      th_1  = th1[6] + th1v[6]*(t - ts[6]) + 0.5*tha_1*pow((t - ts[6]), 2);
    }
  if (ts[7] < t && t < ts[8])
    {
      tha_1 = th1a[8];
      thv_1 = th1v[7] + tha_1*(t - ts[7]);
      th_1  = th1[7] + th1v[7]*(t - ts[7]) + 0.5*tha_1*pow((t - ts[7]), 2);
    }
  if (ts[8] <= t)
    {
    }
  if (4 == step)
    return;
}

inline void
settei(int step)
{
  th1v[1] = 2.0*(th1[1] - th1_i)/t1;
  th1a[1] = th1v[1]/t1;
  th1v[2] = 2.0*(th1[2] - th1[1])/t1 - th1v[1];
  th1a[2] = (th1v[2] - th1v[1])/t1;
  th1v[3] = 2.0*(th1[3] - th1[2])/t1 - th1v[2];
  th1a[3] = (th1v[3] - th1v[2])/t1;
  if (step == 1) return;
  
  th1v[4] = 2.0*(th1[4]-th1[3])/t1 - th1v[3];
  th1a[4] = (th1v[4] - th1v[3])/t1;
  if (step == 2) return;

  th1v[5] = 2.0*(th1[5] - th1[4])/t1 - th1v[4];
  th1a[5] = (th1v[5] - th1v[4])/t1;
  if (step == 3) return;

  th1v[6] = 2.0*(th1[6] - th1[5])/t1 - th1v[5];
  th1a[6] = (th1v[6] - th1v[5])/t1;
  th1[8] = th1_f;
  th1v[8] = 0.0;
  th1v[7] = (th1[8] - th1[6])/t1 - (th1v[6]+th1v[8])/2.0;
  th1a[7] = (th1v[7] - th1v[6])/t1;
  th1[7] = th1[6] + th1v[6] * t1 + th1a[7]*t1*t1/2.0;
  th1a[8] = (th1v[8] - th1v[7])/t1;
  if (step == 4) return;
}

void
Hikaku1(int jj)
{
  int na1, na2, na3;//, na4;
  enemin2 = DBL_MAX;
  for(na3=1; na3<=NN; na3++){
    for(na2=1; na2<=NN; na2++){
      enemin1 = DBL_MAX;
      /* -------------------- */
      for(na1=1; na1<=NN; na1++){
	if (enemin1 > Energ[na1][na2][na3]){
	  enemin1 = Energ[na1][na2][na3];
	  Mat[jj][na2][na3] = na1;
	}
	if (enemin2 > Energ[na1][na2][na3]){
	  enemin2 = Energ[na1][na2][na3];
	  sat[6] = na3;
	  sat[5] = na2;
	  sat[4] = na1;
	}
      }
    }
  }

  //if (DBL_MAX == enemin1)
  //  {
  //    printf("nnn=%d, za=%lf, zv=%lf\n", nnn, r3.za, r3.zv );
  //    printf("enemin1=%lf, DBL_MAX=%lf\n", enemin1, DBL_MAX);
  //    assert(false);
  //  }
}

/*
 * Research the minimam energy
 */
void
tansaku()
{
  reinitialize();
  for(int cur = 1; cur <= 4; cur++){//Magick
    // searching space th1[] 
    for(at[cur+2] = 1; at[cur+2] <= NN; at[cur+2]++){
      for(at[cur+1]=1; at[cur+1] <= NN; at[cur+1]++){
	for(at[cur]=1; at[cur  ] <= NN; at[cur  ]++){
	  
	  th1[cur+2] = Sth1[cur+2] + kiza1*(at[cur+2] - N3);	  
	  th1[cur+1] = Sth1[cur+1] + kiza1*(at[cur+1] - N3);      	  
	  th1[cur  ] = Sth1[cur  ] + kiza1*(at[cur  ] - N3);
	  for(int i = cur; --i > 0;){
	    at[i]  = Mat[i][ at[i+1] ][ at[i+2] ];
	    th1[i] = Sth1[i] + kiza1*(at[i] - N3);
	  }
	  settei(cur);// set array th1v[], th1a[], ntmax;

	  En_sim( ntmax[cur], cur );
	}
      }
    }
  
    Hikaku1(cur);
  }

  return;
}

void
optimal_result(int make_animation_flag)
{
  int cur;
  for(int ii = 6; ii > 3; ii--){ at[ii] = sat[ii]; }//Magick
  for(int ii = 3; ii > 0; ii--){//Magick
     at[ii] = Mat[ii][ at[ii+1] ][ at[ii+2] ];
    sat[ii] = at[ii];
  }  
  for(int ii = 6; ii > 0; ii--){//Magick
    th1[ii] = Sth1[ii] + kiza1*(at[ii] - N3);
  }

  cur = 4;  //Magick
  settei(cur);
  reinitialize();
  ofstream ofs1, ofs2;
  if (make_animation_flag)
    {
      char fname[32];
      sprintf(fname, "dir-dat/saiteki_%05d.dat", nnn);    
      ofs1.open(fname);
      sprintf(fname, "dir-dat/saiteki_ev_%05d.dat", nnn);
      ofs2.open(fname);
    }
  
  //  for(int nt = 1; nt <= ntmax[ cur ]; nt++)
  for(int nt = 0; nt <= ntmax[ cur ]; nt++)  
    {  
      t += dt;
      sim( cur );
      r3.z0 = th_1 ;
      r3.zv = thv_1;
      r3.za = tha_1; 

      En_up(r3);
      //      enemin2 = Ene123 = e3.energy;

      enemin2 = e3.energy;
      
      if (make_animation_flag)
	{
	  graph_standard(ofs1, t,
			 th_1, thv_1, tha_1,
			 e3.ev, e3.ia, e3.tau,
			 e3.energy, e3.ev*e3.ia);
	  
	  graph_volt_constituent(ofs2, t, e3.ev,
				 b1*(r3.zv/r3.radius),
				 b2*(r3.za/r3.radius),
				 b3*(e3.tau));
	  
	  if (10*static_cast<int>(nt*0.1) == nt)
	    {
	      char mfile[32];
	      sprintf(mfile, "dir-dat/dp-%05d.dat", nt);
	      ofstream ofs(mfile);
	      zu(ofs);
	    }
	}// fi make_animation_flag
    } //endfor

  ofs1.close();
  ofs2.close();
  return;
}

void
make_gif_file()
{
  struct {
    char *write;
    char *title;
    int col;
  } const data[8] = {
    { "cmp_z0.gif"    , "z0"    , 2},
    { "cmp_zv.gif"    , "zv"    , 3},
    { "cmp_za.gif"    , "za"    , 4},
    { "cmp_ev.gif"    , "ev"    , 5},
    { "cmp_ia.gif"    , "ia"    , 6},
    { "cmp_tau.gif"   , "tau"   , 7},
    { "cmp_energy.gif", "energy", 8},
    { "cmp_ev_ia.gif" , "ev_ia" , 9},
  };

  FILE *gp = popen("gnuplot", "w");
  fprintf(gp, "set size %lf, %lf\n", 0.6, 0.6 );
  fprintf(gp, "set terminal gif\n");
  fprintf(gp, "set grid\n");

  for (uint ii = 0; ii < sizeof(data)/sizeof(*data); ++ii)
    {
      fprintf(gp, "set key off\n");
      fprintf(gp, "set output 'dir-gif/%s'\n", data[ii].write);
      fprintf(gp, "set title '%s'\n", data[ii].title);

      fprintf(gp, "plot ");
      fprintf(gp, " 'dir-dat/sinc.dat' u 1:%d w l lt 1", data[ii].col);
      for(int jj = 1; jj < REPEAT; ++jj)
	if (plot_flags)
	  fprintf(gp, ",'dir-dat/saiteki_%05d.dat' u 1:%d w l lt 2"
		  , jj, data[ii].col);

      fprintf(gp, ",'dir-dat/saiteki_%05d.dat' u 1:%d w l lw 2 lt 3"
	      , REPEAT, data[ii].col);
      fprintf(gp, ",0  w l lw 1 lt 55");      
      fprintf(gp, "\n");
    }
  
  pclose(gp);
}

ofstream ofs_nnn_kiza("dir-dat/nnn_Sth1_kiza.dat");

void
repeat()
{
  int ff1;
  /* 正弦波の代わりに１回目の探索結果の最適経路を中心位置に選ぶ */
  for(int ii=6; ii > 0; --ii){//Magick
    Sth1[ii] = Sth1[ii] + kiza1 * (sat[ii] - N3);
  }

  ofs_nnn_kiza << nnn <<" "<< kiza1 <<" ";
  for(int ii=1; ii <= 6; ii++){
    ofs_nnn_kiza << Sth1[ii] <<" ";
  }
  ofs_nnn_kiza << std::endl;
  
  {// get summation
    double sum = 0.0;
    for(int ii=6; ii > 0; --ii){//Magick
      sum = sum + fabs(sat[ii] - N3);
    }
    ff1 = static_cast<int>(sum);
  }
  
  if (0 == ff1){
    kiza1 = kiza1 * 0.8;
    printf("kiza1=%lf, nnn=%d\n", kiza1, nnn);
  }
}

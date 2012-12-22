#include <cstdio>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>
#include <fenv.h>
#include "GraphicalPrimitive.hpp"
#include "Manipulator.hpp"
#include "SinCurv.hpp"
#include "Energy.hpp"
#include "debug.h"

using namespace std;

#define g	9.8
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
#define NN	7		// center value

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
// ts[offset], offset=1, Ts=whole time, t1=one part of times
double a1, a2, a3, a4;
double kiza1, kiza2, kiza3;

int sat[1 + 6], at[1 + 6];
int nnn; // repeat times

// using DP_function
void repeat ();
void tansaku(int make_animation_flag);
void makeAnimetion();
void settei (int step);
void sim(int step);
double tha_1, thv_1, th_1;	// return sim() parametar
double tha_2, thv_2, th_2;	// return sim() parametar
double tha_3, thv_3, th_3;	// return sim() parametar
void Hikaku1 (int j);
void zu(ofstream &ofs);
void sinc();

void reinitialize()
{
  t = 0.0;
 
  e1 = Energy();
  r1.th = 0.0; r1.thv = 0.0; r1.tha = 0.0;

  e2 = Energy();
  r2.th = 0.0; r2.thv = 0.0; r2.tha = 0.0;
   
  e3 = Energy();
  r3.z0 = 0.0; r3.zv = 0.0; r3.za = 0.0;
}

void init()
{
  // \< Don't change follow lines sequence.
  const struct Point{
    double x, y, z;
  } pnt_i  = { 0.05        , 0.10       , 0.0   } // real end_effector start
  , pnt_f  = {-0.05        , 0.10       , 0.0   } // real end_effector end
  , pnt_v0 = {-0.05 - 0.005, 0.10 + 0.05, 0.0   } // virtual base point
  , pnt_r0 = { 0.00        , 0.00       , 0.0618}; // real base point
  
  const double
    vth_i     = atan2(pnt_i.y - pnt_v0.y, pnt_i.x - pnt_v0.x),
    vth_f     = atan2(pnt_f.y - pnt_v0.y, pnt_f.x - pnt_v0.x),
    rth12     = atan2(pnt_i.y - pnt_r0.y, pnt_i.x - pnt_r0.x),
    line_i_r0 = hypot(pnt_i.x - pnt_r0.x, pnt_i.y - pnt_r0.y),
    line_i_v0 = hypot(pnt_i.x - pnt_v0.x, pnt_i.y - pnt_v0.y);  

  const double rth1  = rth12 - acos(0.5*line_i_r0/L1);

  struct Point pnt_r1 = {
    pnt_r0.x + L1*cos(rth1), pnt_r0.y + L1*sin(rth1), pnt_r0.z,    
  };
  
  //  double rth2  = 2.0*acos(0.5*line_i_r0/L);
  const double rth2 = atan2(pnt_i.y - pnt_r1.y, pnt_i.x - pnt_r1.x);

  const struct Config {
    double ini, fin;// using th1[], th2[], th3[]
  } cnf_th1 = { 0.0   , 0.0618 } // distance
  , cnf_th2 = { vth_i , vth_f  } // degree
  , cnf_th3 = { 0.0618, 0.0    }; // distance
  // -- Don't change sequence \>
  
  const double radius = 0.04;	// pulley's radius
  // optimal pulley's radius = `r = 0.04` when rouund-trip.
  
  // initialize global value
  t  = 0.0; Ts = 0.64; t1 = Ts/8.0;// pull pulley_test2.cpp 
  //  t1 = Ts/(double)(DIV)
  //  for(int i=0; i <= DIV; i++){ ts[i] = t1*i; }
  for(int i=0; i <= 8; i++){ ts[i] = t1*i; }
  
  b1 = Kv + Ra*Dm/Kt; b2 = Ra*Im/Kt; b3 = Ra/Kt;
  
  //  N3 = static_cast<int>(NN * 0.5 + 0.5);
  N3 = static_cast<int>(NN * 0.5 + 0.5);

  //kiza1 = 0.001;
  kiza1 = 0.0001;
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
  //  sc1 = SinCurv(0.64,   0.0, 0.0618, 0.0, 0.0);// 上昇 pull lab/work2.cpp
  //  sc3 = SinCurv(0.64,   1.0,    0.0, 0.0, 0.0);// 下降
  
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
    ofstream ofs("sin_th1.dat");
    for(int i=0; i <= 8; ++i){
      Sth1[i] = sc1.sin_th1(ts[i]);
      ofs << ts[i] <<" "<< Sth1[i] << std::endl;
    }
    ofs.close();
  }
  
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

void view_animation( std::string prefix )
{
  FILE *gp = popen("gnuplot","w");
  fprintf(gp, "set xrange[%lf:%lf]\n", -0.08, 0.10);
  fprintf(gp, "set yrange[%lf:%lf]\n", -0.02, 0.16);
  fprintf(gp, "set zrange[%lf:%lf]\n",  0.00, 0.20);
  fprintf(gp, "set view %lf, %lf\n", 67.0, 133.0);
  for(int ii=10; ii < static_cast<int>(ts[8]/dt); ii+=10){
    fprintf(gp, "splot 'dir-dat/%s-%05d.dat' w l\n"
	    , prefix.c_str(), ii);
    fprintf(gp, "pause 0.01\n");      
  }
  pclose(gp);
}

void zu()
{
  cout <<__FILE__<<"_"<<__FUNCTION__<<"_"<<__LINE__<< endl;
}

void zu(std::ofstream& ofs)
{
  cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<< endl;
  Circle3D real1_p0(r1.x0, r1.y0, r1.z0, 0.01);
  Circle3D real1_p1(r1.x , r1.y , r1.z0, 0.01);

  Circle3D real2_p0(r2.x0, r2.y0, r2.z0, 0.01);
  Circle3D real2_p1(r2.x , r2.y , r2.z0, 0.01);
  Circle3D real2_p2(r3.x0, r3.y0, r3.z0, 0.01);

  Circle3D real3_p0(r3.x0, r3.y0, r3.z0, 0.01);
  Circle3D real3_p1(r3.x0, r3.y0, r3.z , 0.01);
  
  Line3D real1_ln(r1.x0, r1.y0, r1.z0, r1.x, r1.y, r1.z0 );
  Line3D real2_ln(r2.x0, r2.y0, r1.z0, r2.x, r2.y, r1.z0 );
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

void En_up(Prismatic& r3)
{
  double F = r3.m*(r3.za + g);
  e3.tau   = F*r3.radius;  	// pulley's radius and torque.
  e3.ev    = b1*(r3.zv/r3.radius) + b2*(r3.za/r3.radius) + b3*e3.tau;
  e3.ia    = (e3.ev - Kv * r3.zv/r3.radius)/Ra;

  if (e3.ev*e3.ia > 0.0){
    e3.energy = e3.energy + e3.ev*e3.ia*dt;
    //Nov5 11:20    printf("line=%d, e3.energy=%lf\n", __LINE__, e3.energy);
  }
}

void En_sim(int ntmax, int cur)
{
  reinitialize();
  for(int nt = 1; nt <= ntmax; nt++){
    t += dt;
    sim(cur);// set th_1, thv_1, tha_1
    r3.z0 = th_1 ; r3.zv = thv_1; r3.za = tha_1; 
    En_up(r3);
  }
  Energ[at[cur]][at[cur+1]][at[cur+2]] = Ene123 = e3.energy;
}

void En_direct(Revolute& r1, Revolute& r2)
{
  std::cout << "Nothing do this function. " << __FUNCTION__ << std::endl;
}

void En_down(Prismatic& r3)
{
  std::cout << "Nothing do this function. " << __FUNCTION__ << std::endl;
  double F = r3.m*(r3.za + g);
  e1.tau= F*r3.radius;  	// pulley's radius and torque.
  e1.ev = b1*(r3.zv/r3.radius) + b2*(r3.za/r3.radius) + b3*e1.tau;
  e1.ia = (e1.ev - Kv * r3.zv/r3.radius)/Ra;
  if (e1.ev*e1.ia > 0.0){
    e1.energy = e1.energy + e1.ev*e1.ia*dt;
  }  
}

int
main(int argc, char *argv[])
{
  feenableexcept (FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  Manipulator m(m1, LL, 0.0, 0.0, 0.0);
  // Manipulator(mass, length, x0, y0, z0);
  init();
  sinc();
  view_animation("e1");		// view sinc animation
  
  {// plot enemin2 block
    ofstream ofs("dir-dat/enemin2.dat");
    //for (int nnn = 1; nnn <= 13; nnn++) {
    for (nnn = 1; nnn <= 60; nnn++){
      repeat();
      enemin2 = 0.0;

      tansaku(0);

      Enemin2[nnn] = enemin2;

      if ( (nnn > 15) && (Enemin2[nnn] > Enemin2[nnn - 1]) )
	kiza1 *= 1.2;
      
      ofs << nnn <<" "<< enemin2 <<" "<< kiza1 << endl;
    }
  }

  repeat ();
  enemin2 = 0.0;
  tansaku(1);
  // if speed is more precious replace this function `make_animation`.
  
  view_animation("dp");
  
  return 0;
}



void sinc()
{
  ofstream ofs;
  ofs.open("dir-dat/sinc.dat");

  reinitialize();
  for (int i = 0; i < static_cast <int>(ts[8]/dt); ++i){
    t += dt;
    th_1  = sc1.sin_th1(t);
    thv_1 = sc1.sin_th1v(t);
    tha_1 = sc1.sin_th1a(t);

    r3.z0 = th_1 ; r3.zv = thv_1; r3.za = tha_1;
    
    En_up(r3);

    if (i % 10 == 0){
      char fname[32];
      sprintf(fname, "dir-dat/e1-%05d.dat", i);
      ofstream fout(fname);
      zu(fout);
    }

    ofs << t         << setw(16)
	<< th_1      << setw(16)
	<< thv_1     << setw(16)
	<< tha_1     << setw(16)
	<< e3.ev     << setw(16)
	<< e3.ia     << setw(16)
	<< e3.tau    << setw(16)
	<< e3.energy << setw(16)            
	<< std::endl;
  }
  ofs.close();
}

void
sim(int step)
{
  if (t < ts[1]) {
    tha_1 = th1a[1];
    thv_1 = tha_1 * t;
    th_1  = 0.5 * tha_1 * t * t + th1_i;
  }
  if (ts[1] < t && t < ts[2]) {
    tha_1 = th1a[2];
    thv_1 = th1v[1] + tha_1 * (t - ts[1]);
    th_1  = th1[1] + th1v[1]*(t - ts[1]) + 0.5*tha_1*pow((t - ts[1]), 2);
  }
  if (ts[2] < t && t < ts[3]){
    tha_1 = th1a[3];
    thv_1 = th1v[2] + tha_1 * (t - ts[2]);
    th_1  = th1[2] + th1v[2]*(t - ts[2]) + 0.5*tha_1*pow((t - ts[2]), 2);
  }

  if (step == 1)
    return;

  if (ts[3] < t && t < ts[4]){
    tha_1 = th1a[4];
    thv_1 = th1v[3] + tha_1*(t - ts[3]);
    th_1 = th1[3] + th1v[3]*(t - ts[3]) + 0.5 * tha_1 * pow ((t - ts[3]), 2);
  }
  if (step == 2)
    return;

  if (ts[4] < t && t < ts[5]) {
    tha_1 = th1a[5];
    thv_1 = th1v[4] + tha_1 * (t - ts[4]);
    th_1  = th1[4] + th1v[4] * (t - ts[4])
      + 0.5 * tha_1 * pow ((t - ts[4]), 2);
  }
  if (step == 3)
    return;

  if (ts[5] < t && t < ts[6]) {
    tha_1 = th1a[6];
    thv_1 = th1v[5] + tha_1 * (t - ts[5]);
    th_1 = th1[5] + th1v[5]*(t - ts[5]) + 0.5*tha_1*pow((t - ts[5]), 2);
  }
  if (ts[6] < t && t < ts[7]) {
    tha_1 = th1a[7];
    thv_1 = th1v[6] + tha_1 * (t - ts[6]);
    th_1 = th1[6] + th1v[6]*(t - ts[6]) + 0.5*tha_1*pow((t - ts[6]), 2);
  }
  if (ts[7] < t && t < ts[8]) {
    tha_1 = th1a[8];
    thv_1 = th1v[7] + tha_1*(t - ts[7]);
    th_1 = th1[7] + th1v[7]*(t - ts[7]) + 0.5*tha_1*pow((t - ts[7]), 2);
  }
  if (ts[8] <= t) {
  }
  if (step == 4)
    return;
}

void
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
  enemin2 = 1000000.0;
  for(na3=1; na3<=NN; na3++){
    for(na2=1; na2<=NN; na2++){
      enemin1 = 1000000.0;
      /* -------------------- */
      for(na1=1; na1<=NN; na1++){
	if (enemin1 > Energ[na1][na2][na3]){
	  enemin1 = Energ[na1][na2][na3];
	  Mat[jj][na2][na3] = na1;
	}
	if (enemin2 > Energ[na1][na2][na3]){
	  enemin2 = Energ[na1][na2][na3];
	  sat[6] = na3; sat[5] = na2; sat[4] = na1;
	}
      }
    }
  }
}

void
tansaku(int make_animation_f)
{
  int cur;
  int step[1+4] = {// nt's max length
    0,
    static_cast<int>(ts[3]/dt),
    static_cast<int>(ts[4]/dt),
    static_cast<int>(ts[5]/dt),
    static_cast<int>(ts[8]/dt),    
  };
  //for(uint i=0; i < sizeof(step)/sizeof(*step); ++i){
  //  step[i] = static_cast<int>((ts[i+2])/dt);
  //}
  //step[sizeof(step)/sizeof(*step) - 1]
  //  = static_cast<int>( (ts[8])/dt );

  reinitialize();
  for(cur = 1; cur <= 4; cur++){
    // searching space th1[] 
    for(at[cur+2] = 1; at[cur+2] <= NN; at[cur+2]++){
      for(at[cur+1]=1; at[cur+1] <= NN; at[cur+1]++){
	for(at[cur]=1; at[cur  ] <= NN; at[cur  ]++){
	  
	  th1[cur+2] = Sth1[cur+2] + kiza1*(at[cur+2] - N3);	  
	  th1[cur+1] = Sth1[cur+1] + kiza1*(at[cur+1] - N3);      	  
	  th1[cur  ] = Sth1[cur  ] + kiza1*(at[cur  ] - N3);
	  //--------------------------------------------------
	  for(int i = cur; --i > 0;){
	    at[i]  = Mat[i][ at[i+1] ][ at[i+2] ];
	    th1[i] = Sth1[i] + kiza1*(at[i] - N3);
	  }
	  // --------------------------------------------------

	  settei(cur);// set array th1v[], th1a[], step;
	  En_sim( step[cur], cur );
	}
      }
    }
    
    Energ[ at[cur] ][ at[cur+1] ][ at[cur+2] ] = Ene123;

    Hikaku1(cur);
  }

  for(int ii=6; ii > 3; ii--){ at[ii] = sat[ii]; }
  for(int ii = 3; ii > 0; ii--){
    at[ii] = Mat[ii][ at[ii+1] ][ at[ii+2] ];
    sat[ii] = at[ii];
  }  
  for(int ii=6; ii > 0; ii--){
    th1[ii] = Sth1[ii] + kiza1*(at[ii] - N3);
  }

  // simulate optimal route's
  cur = 4;			// if divided number changed
  settei(cur);
  reinitialize();
  {// draw "dir-dat/saiteki_%05d.dat"
    char fname[32];
    sprintf(fname, "dir-dat/saiteki_%05d.dat", nnn);    
    ofstream ofs1("t_th.dat"), ofs2("t_e3.dat"), ofs3(fname);
  
    for(int nt = 1; nt <= step[ cur ]; nt++){  
      t += dt;
      sim( cur );
      r3.z0 = th_1 ; r3.zv = thv_1; r3.za = tha_1; 
      En_up(r3);
      
      enemin2 = Ene123 = e3.energy;
    
      ofs1 << t      <<" "<< th_1   <<" "
	   << thv_1  <<" "<< tha_1  <<" "
	   << Ene123 <<" "<< std::endl;

      ofs2 << t         <<" "<< e3.ev     <<" "
	   << e3.ia     <<" "<< e3.tau    <<" "
	   << e3.energy <<" "<< std::endl;

      ofs3 << t <<" "
	   << e3.energy <<" "
	   << th_1 <<" "
	   << thv_1 <<" "
	   << tha_1 << std::endl;
      
      if (make_animation_f){
	if (10*static_cast<int>(nt*0.1) == nt){
	  char mfile[32];
	  sprintf(mfile, "dir-dat/dp-%05d.dat", nt);
	  ofstream ofs(mfile);
	  zu(ofs);
	}
      }
    }
  } // draw "dir-dat/saiteki_%05d.dat

  return;
}

void makeAnimetion()
{
  cout << __FILE__     <<"_"<< __FUNCTION__ <<"_"<< __LINE__ <<endl;
}

void repeat()
{
  int ff1;
  printf("nnn=%3d\n", nnn);
  /* 正弦波の代わりに１回目の探索結果の最適経路を中心位置に選ぶ */
  for(int ii=6; ii > 0; --ii){
    Sth1[ii] = Sth1[ii] + kiza1 * (sat[ii] - N3);
  }

  {// get summation
    double sum = 0.0;
    for(int ii=6; ii > 0; --ii){
      sum = sum + fabs(sat[ii] - N3);
    }
    ff1 = static_cast<int>(sum);
  }
  
  if (ff1 == 0){
    kiza1 = kiza1 * 0.8;
    printf("kiza1=%lf, nnn=%d\n", kiza1, nnn);
  }
}

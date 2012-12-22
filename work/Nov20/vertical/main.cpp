#include <cstdio>
#include <cmath>
#include <cerrno>
#include <cfloat>
#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <sys/times.h>
#include "config.hpp"
#include "Manipulator.hpp"
//------------------ advanced course of robotics -------------
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
#define Lg1	 0.04
#define Lg2	 0.04
#define Lg3	 0.04
#define IG1	 1.73e-05   
#define IG2	1.73e-05
#define m1	0.0202
#define m2	0.0202
#define m3	0.0202
#define NN	7		// center value
//#define NN	3		// center value
#define DIV	16		// division number

using namespace std;

struct Return_sim_value_t th_;
//struct Manipulator v1, v2, r1, r2, r3;

double b1, b2, b3;		// using moter
double a1, a2, a3, a4, a5;	// inerita?
double t, Ts, t1, ts[20], ttf, tt, Pdeg;
double th1_i, th1_f, th2_i, th2_f, fai_2i, fai_2f;
double th1v_i;
double ac1, Bc1, mL;
double omg1, Amp1, omg2, Amp2;	// using sincurve

//namespace energy{
double ev1, ev2, ev3, ia1, ia2, ia3;
double tau_1, tau_2, tau_3, fz_3;// using En1, En2, En3
double Ene_1, Ene_2, Ene_3, enemin1, enemin2;
double Ene_123;
//} 
// virtual manipulator horizontal
//struct Return_sim_value_t {double s, v, a;} th_;

double Sth1[20], kiz1a, th_1, thv_1, tha_1, th1[20], th1v[20], th1a[20];
double Sth2[20], kiz2a, th_2, thv_2, tha_2, th2[20], th2v[20], th2a[20];
// virtual manipulator vertical
double Energ[15][15][15], Enemin2[50];
//double cth1, cth2, cth12, sth1, sth2, sth12; local register
//double detJ, XXa, YYa, PL3;

double fxx,fyy,fxa,fya,keirosokudo,keirosokudo2;
//int n, nnn, i;
int nt, nt11, nt22, nt33;
int Mat[13][13][13], N3;
int sat[1+DIV], at[1+DIV];
int flag11, ff1, ff2, flag34;

struct Environment{
  Conf_value *ini, *fin;
  Linear_curve_t *fk1, *fk2, *fkv;
  Sin_curve *sc1, *scv;
  BallScrew *bs;
  Manipulator *v1, *v2, *r1, *r2, *r3;
  Return_sim_value_t *th_;
};

void sinc(Environment &env);
double En1(const Return_sim_value_t &th_);
//double En_vertical(const Return_sim_value_t &th_);
double En_prism(const struct Return_sim_value_t &th_,
		Linear_curve_t& fk1, Linear_curve_t& fk2,
		Conf_value& ini, Conf_value& fin,		
		Manipulator& v1,
		Manipulator& r1, Manipulator& r2);
double En_vertical(const struct Return_sim_value_t &z_,
		   BallScrew& bs, 
		   Manipulator &v1, Manipulator &r3);

void zu(FILE *_fp1, FILE *_fp2, FILE *_fp3, Environment &e);
void zu_virtual(FILE *_fp, Manipulator& v1, Return_sim_value_t& th_);
void zu_real(FILE *_fp,
	     Manipulator& r1,Manipulator& r2,Manipulator& r3,
	     Return_sim_value_t& th_);
void zu_pathway(FILE *_fp, Manipulator& v1);
  
void draw_circle();
void draw_circle_fp(FILE *_fp, double x, double y, double z, double r);

inline double det(double a11, double a12, double a21, double a22)
{
  return a11*a22 - a12*a21;
}
//**********************************************************************
void init_variable(FILE *_fp, Environment &env)
{
  FILE *fp1 = fopen("VirtualManipulator.dat","w");
  FILE *fp2 = fopen("RealManipulator.dat","w");
  double l_hypot, l_base;
  struct Sin_curve *psc1, *pscv;
  class BallScrew *pbs;
  struct Manipulator
    *pv1 = new Manipulator(0.0, L1), *pv2 = new Manipulator(0.0, 0.0),
    *pr1 = new Manipulator(m1 , L1), *pr2 = new Manipulator(m2 ,  L2),
    *pr3 = new Manipulator(m3 , L3);
  pr3->Ig = 0.0;		// !!exception case!!
  struct Conf_value
    *pini = new Conf_value(),*pfin = new Conf_value();
  struct Linear_curve_t *pfk1, *pfk2, *pfkv;
  struct Manipulator &v1=*pv1, &v2=*pv2, &r1=*pr1, &r2=*pr2, &r3=*pr3;
  //-initialize value
  // setting time variable
  Ts = 0.64 ;
  t1 = Ts/static_cast<double>(DIV); //t1 = Ts/16.0;

  for(int i = 1; i <= DIV; i++){
    ts[i] = t1*i; fprintf(_fp,"%d ts[%d]=%lf\n", __LINE__, i, ts[i]);
  }
  ttf = ts[DIV];
  
  // actuator
  b1 = Kv + Ra*Dm/Kt; b2 = Ra*Im/Kt; b3 = Ra/Kt; // dt=0.001;
  v1.setP0( -0.05 - 0.005,  0.10 + 0.05, 0.0);
  {// start and end point
    pini->x  =  0.05; pini->y  = 0.10; pini->z  =  0.00;
    pfin->x  = -0.05; pfin->y  = 0.10; pfin->z  =  0.0618;

    l_hypot = hypot(pini->x - v1.x0, pini->y  - v1.y0);
    l_base  = hypot(pfin->x - v1.x0, pfin->y - v1.y0);

    pini->vth = atan2(pini->y - v1.y0, pini->x - v1.x0);
    pfin->vth = atan2(pfin->y - v1.y0, pfin->x - v1.x0);

    l_hypot = hypot(pini->x, pini->y);
    pini->rth1 = atan2(pini->y, pini->x) - acos(0.5*l_hypot/r1.l);
    pini->rth2 = 2.0*acos(0.5*l_hypot/r1.l);
  }// start and end point
  th_.s = pini->vth;
  // virtual manipulator
  v1.m = m1;
  v1.th = th_.s;
  //  v1.l = 0.05; v1.r = v1.l/2.0;
  v1.x1 = pini->x;
  v1.y1 = pini->y;

  
  v2 = v1;

  //r1.x0 = 0.0;  r1.y0 = 0.0;  r1.z0 = HH;
  r1.setP0(0.0, 0.0, HH);	// x, y, z
  r1.x  = r1.x0 + r1.r*cos(pini->rth1);
  r1.x1 = r1.x0 + r1.l*cos(pini->rth1);
  r1.y  = r1.y0 + r1.r*sin(pini->rth1);
  r1.y1 = r1.y0 + r1.l*sin(pini->rth1);
  r1.z  = HH;
  r1.z1 = HH;
  r1.th = pini->rth1;  
  //  r2.x0 = r1.x1;  r2.y0 = r1.y1;  r2.z0 = HH;
  r2.setP0(r1.x1, r1.y1, HH);
  r2.x  = r2.x0 + r2.r*cos(pini->rth1 + pini->rth2);
  r2.x1 = r2.x0 + r2.l*cos(pini->rth1 + pini->rth2);
  r2.y  = r2.y0 + r2.r*sin(pini->rth1 + pini->rth2);
  r2.y1 = r2.y0 + r2.l*sin(pini->rth1 + pini->rth2);  
  r2.z = HH;
  r2.z1 = HH;
  r2.th = pini->rth2;
  //  r3.x0 = r2.x1; r3.y0 = r2.y1;  r3.z0 = HH   ;
  r3.setP0(r2.x1, r2.y1, HH);
  r3.x  = r2.x1; r3.x1 = r2.x1;
  r3.y  = r2.y1; r3.y1 = r2.y1;
  r3.z  = r3.z0 - r3.l/2.0;
  r3.z1 = r3.z0 - r3.l;
  r3.th = 0.0;

  //  fk1.a = 0.0; fk1.b = 0.1;
  //  fk2.a = tan(th_.s); fk2.b = v1.y0 - fk2.a*v1.x0;
  pfk1 = new Linear_curve_t(0.0, 0.1);
  pfk2 = new Linear_curve_t(tan(th_.s), v1.y0 - tan(th_.s)*v1.x0);

  N3 = static_cast<int>(NN*0.5);
  //------------- Scara Robot --------------------------
  a1 = r1.Ig + r1.m*r1.r*r1.r + r2.m*r1.l*r1.l + r3.m*r1.l*r1.l;
  a2 = r2.Ig + r2.m*r2.r*r2.r + r3.m*r2.l*r2.l;
  a3 = r3.Ig;
  a4 = r2.m*r1.l*r2.r + r3.m*r1.l*r2.l;
  //a5 = ;
  //-----------------------------------------------------
  {// sono uti shusei
    th1_i = pini->vth; th1_f = pfin->vth;
    th1v_i = 0.0 ;
    th1[0]   = th1_i ; th1v[0]   = th1v_i;
    th1[DIV] = th1_f ; th1v[DIV] = 0.0;
  }
  // libration and amplitude
  Pdeg = 180.0/M_PI ;		// may be ballscrew. bad he did not reach.
  omg1 = M_PI/ttf;
  Amp1 = ((th1_f - th1_i)/2.0)*omg1;
  //Amp1 = ((fin.vth - ini.vth)/2.0)*omg1;
  ac1 = (Amp1 - th1v_i)/ttf;
  Bc1 = (M_PI/2.0)*((th1_f - th1_i)/ttf - (Amp1 - th1v_i)/2.0);

  psc1 = new Sin_curve( ttf, pini->vth, pfin->vth, 0.0, 0.0); 
  // ttf(:time to final), (th1), (th1v)
  pscv = new Sin_curve( ttf, pini->z  , pfin->z  , 0.0, 0.0);
  pbs = new BallScrew(Im, r1.m, 2.0*M_PI);
 // motor side: moment of inerita 
  for(int n = 1; n <= DIV; n++){
    Sth1[n] = psc1->sin_th1(ts[n]);
    printf("%s-%d\n", __FUNCTION__, __LINE__);
    psc1->print_value(stdout, ts[n]);
    fprintf(_fp,"%d tt=%lf,Sth1[%d]=%lf %lf\n"
	    ,__LINE__, tt, n, Sth1[n], ts[n]);
  }

  kiz1a = 0.04; kiz2a = 0.04;

  //  zu_virtual(fp1);
  //  zu_real(fp2);
  *pv1 = v1; *pv2 = v2;
  *pr1 = r1; *pr2 = r2; *pr3 = r3;
  cout <<__FUNCTION__<<"-"<<__LINE__ << endl;
  {// setup environment
    env.ini = pini; env.fin = pfin;
    env.fk1 = pfk1; env.fk2 = pfk2; env.fkv = pfkv;
    env.sc1 = psc1; env.scv = pscv;
    env.bs  = pbs;
    env.v1 = &v1; env.v2 = &v2;
    env.r1 = &r1; env.r2 = &r2; env.r3 = &r3;
    env.th_ = &th_;
  } // setup environment
  
  {//fprintf_line
    fprintf(_fp,"%d Ts=%lf,t1=%lf,DIV=%d\n",__LINE__, Ts, t1, DIV);  
    fprintf(_fp,"%d b1=%lf,b2=%lf,b3=%lf\n",__LINE__,b1,b2,b3);
    fprintf(_fp, "%d v1.x0=%lf,v1.y0=%lf,v1.z0=%lf\n"
	    ,__LINE__, v1.x0, v1.y0, v1.z0);
    fprintf(_fp,"%d ini.vtht=%lf,(d)%lf\n"
	    ,__LINE__, pini->vth, pini->vth*180.0/M_PI);
    fprintf(stdout,"%d ini.vtht=%lf,(d)%lf\n"
	    ,__LINE__, pini->vth, pini->vth*180.0/M_PI);
    fprintf(_fp,"%d fin.vth=%lf,(d)%lf\n"
	    ,__LINE__, pfin->vth, pfin->vth*180.0/M_PI);
    fprintf(stdout,"%d fin.vth=%lf,(d)%lf\n"
	    ,__LINE__, pfin->vth, pfin->vth*180.0/M_PI);        
    fprintf(_fp,"%d hypot(%lf,%lf)=%lf,/2=%lf\n"
	    ,__LINE__, pini->x, pini->y, l_hypot, l_hypot/2.0);    
    fprintf(_fp,"%d rth1.init=%lf,(d)%lf\n"
	    ,__LINE__, pini->rth1, pini->rth1*180.0/M_PI);    
    fprintf(_fp,"%d rth2.init=%lf,(d)%lf\n"
	    ,__LINE__, pini->rth2, pini->rth2*180.0/M_PI);  
    fprintf(_fp,"%d v1.m=%lf,v1.th=%lf\n",__LINE__, v1.m, v1.th);  
    fprintf(_fp, "%d v1.l=%lf,v1.r=%lf\n",__LINE__, v1.l, v1.r);
  
    fprintf(_fp,"%d fk1.a=%lf,fk1.b=%lf\n", __LINE__, pfk1->a, pfk1->b);
    fprintf(_fp,"%d fk2.a=%lf,fk2.b=%lf\n", __LINE__, pfk2->a, pfk2->b);
    fprintf(_fp,"%d th_.s=%8.4f tan(th_.s)=%8.4f\n"
	    , __LINE__, th_.s, tan(th_.s));
    fprintf(_fp,"%d cfg_rthi=%lf,fth2_i=%lf\n"
	    , __LINE__, pini->rth1, pini->rth2);  
    fprintf(_fp,"%d NN=%3d  N3=%3d\n",__LINE__,NN,N3);
    fprintf(_fp,"%d hypot=%lf,base=%lf,l_base/l_hypot=%lf\n"
	    ,__LINE__, l_hypot, l_base, l_base/l_hypot);
    v1.print_value(_fp,"v1"); r1.print_value(_fp,"r1");
    r2.print_value(_fp,"r2"); r3.print_value(_fp,"r3");  
    fprintf(_fp,"%d r1.Ig=%lf,r1.l=%lf,r1.r=%lf\n"
	    ,__LINE__, r1.Ig, r1.l, r1.r);
    fprintf(_fp,"%d r2.Ig=%lf,r2.l=%lf,r2.r=%lf\n"
	    ,__LINE__, r2.Ig, r2.l, r2.r);
    fprintf(_fp,"%d r3.Ig=%lf,r3.l=%lf,r3.r=%lf\n"
	    ,__LINE__, r3.Ig, r3.l, r3.r);
    fprintf(_fp,"%d a1=%lf\n",__LINE__, a1);
    fprintf(_fp,"%d a2=%lf\n",__LINE__, a2);
    fprintf(_fp,"%d a3=%lf\n",__LINE__, a3);
    fprintf(_fp,"%d a4=%lf\n",__LINE__, a4);
    fprintf(_fp,"%d IG1=%lf,IG2=%lf,IG3=%lf\n"
	    ,__LINE__, r1.Ig, r2.Ig, r3.Ig);
    fprintf(_fp,"%d m1=%lf,m2=%lf,m3=%lf\n",__LINE__, r1.m, r2.m, r3.m);
    fprintf(_fp,"%d L1=%lf,L2=%lf,L3=%lf\n",__LINE__, r1.l, r2.l, r3.l);
    fprintf(_fp,"%d Lg1=%lf,Lg2=%lf,Lg3=%lf\n"
	    ,__LINE__, r1.r, r2.r, r3.r);    
    fprintf(_fp,"%d th1_i=%lf, th1_f=%lf\n",__LINE__, th1_i, th1_f);
    fprintf(_fp,"%d th1_i(d)=%lf, th1_f(d)=%lf\n",__LINE__
	    , th1_i*180.0/M_PI, th1_f*180.0/M_PI);  
    fprintf(_fp,"%d Amp1=%lf, th1v_i=%lf\n",__LINE__, Amp1, th1v_i);
    fprintf(_fp,"%d Pdeg=%lf,omg1=%lf\n",__LINE__, Pdeg, omg1);  
    fprintf(_fp,"%d ac1=%lf,ttf=%lf,Bc1=%lf\n",__LINE__,ac1,ttf,Bc1);
    fprintf(_fp,"%d kiz1a=%lf\n",__LINE__,kiz1a);
    fprintf(_fp,"%d kiz2a=%lf\n",__LINE__,kiz2a);
    fprintf(_fp,"%d nt22=%d\n",__LINE__, nt22);
    fprintf(_fp,"%d nt11=%d\n",__LINE__, nt11);    
  }//fprintf_line
  
  fclose(fp1); fclose(fp2);
}

double gettimeofday_sec()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int main(void)
{
  double start, end;
  Environment env;
  FILE *fp1;
  fp1 = fopen("work2_variable.dat","w");

  start = gettimeofday_sec();

  cout << __FUNCTION__<<"-"<<__LINE__ << endl;
  
  init_variable(fp1, env);
  cout << __FUNCTION__<<"-"<<__LINE__ << endl;
  
  fclose(fp1);
  sinc(env);
  end = gettimeofday_sec();
  //Animation();
  cout << __LINE__ << endl;
  
  cout << "&env.v1" << (env.v1) << endl;

  cout << __LINE__ << endl;
  printf("time=%g,start=>%g,end=%g\n", end - start, start, end);
  return 0 ;
}
//************************************************************
void sinc(Environment& env)
{
  char mfile[32];
  FILE *gp1, *gp2;
  Manipulator& r1 = *(env.r1);
  Manipulator& r2 = *(env.r2);
  Manipulator& r3 = *(env.r3);

  Sin_curve& scv = *(env.scv);
  Sin_curve& sc1 = *(env.sc1);
  
  ofstream ofs0("dir-dat/sinc1.dat");
  ofstream ofs1("dir-dat/end_effector.dat");
  ofstream ofs4("dir-dat/singular.dat");
  ofstream ofs8("dir-dat/virtual_joint.dat");
  ofstream ofs9("dir-dat/fxx_fyy.dat");
  ofstream ofs_ene1("dir-dat/energy1.dat");
  ofstream ofs_ene2("dir-dat/energy2.dat");
  ofstream ofs_ene3("dir-dat/energy3.dat");

  gp1 = fopen("dir-dat/trjc_real.dat","w");
  gp2 = fopen("dir-dat/trjc_virtual.dat","w");
  
  t = 0.0 ; Ene_1 = 0.0 ; Ene_2 = 0.0; 
  //  r1.th  = ini.rth1;  r1.thv = 0.0;
  //  r2.th  = ini.rth2;  r2.thv = 0.0;

  for(nt = 0; nt <= static_cast<int>(ttf/dt); nt++){
    t = t + dt;       
//    th_.s = sc1->sin_th1(t);
//    th_.v = sc1->sin_th1v(t);
//    th_.a = sc1->sin_th1a(t);

    // this value give [m]
    th_.s = scv.sin_th1(t);
    th_.v = scv.sin_th1v(t);
    th_.a = scv.sin_th1a(t);
    //En1(th_);
    //En_vertical(th_);
    //En_vertical(th_);
    En_vertical(*(env.th_),	// latter this value replacement z_
		*(env.bs), 
		*(env.v1), *(env.r3));
    //    cout << env.th_ <<" "<< &th_ << endl;
    //    assert(false);
//    En_prism(*(env.th_),
//	     *(env.fk1), *(env.fk2),
//	     *(env.ini), *(env.fin),
//	     *(env.v1),
//	     *(env.r1), *(env.r2));
    //ofstream ofs0("dir-dat/sinc1.dat");    
    ofs0 << setw(16) << t      
	 << setw(16) << r1.th  << setw(16) << r1.thv << setw(16) << r1.tha
	 << setw(16) << r2.th  << setw(16) << r2.thv << setw(16) << r2.tha
	 << setw(16) << r3.th  << setw(16) << r3.thv << setw(16) << r3.tha
	 << setw(16) << r3.z1  << setw(16) << r3.zv  << setw(16) << r3.za
	 << endl;
    //  ofstream ofs_ene1("dir-dat/energy1.dat");
    ofs_ene1 << t   << setw(16)  << Ene_1 << setw(16)
	     << ev1 << setw(16)  << ia1   << setw(16)
	     << tau_1 << Ene_123 << endl;
    //  ofstream ofs_ene2("dir-dat/energy2.dat");
    ofs_ene2 << t   << setw(16)  << Ene_2 << setw(16)
	     << ev2 << setw(16)  << ia2   << setw(16)
	     << tau_2 << Ene_123 << endl;
    //  ofstream ofs_ene3("dir-dat/energy3.dat");
    ofs_ene3 << t   << setw(16)  << Ene_3 << setw(16)
	     << ev3 << setw(16)  << ia3   << setw(16)
	     << tau_3 << Ene_123 << endl;
    //ofstream ofs4("dir-dat/singular.dat");
    ofs4 << t << setw(16)
	 << tan(th_.s) << setw(16)
	 << th_.s*180.0/M_PI << setw(16)
	 << th_.v << setw(16)
	 << th_.a << endl;
    //  ofstream ofs8("dir-dat/virtual_joint.dat");
    ofs8 << t << setw(16)
	 << th_.s << setw(16)
	 << th_.v << setw(16)
	 << th_.a << setw(16)
	 << th_.s*180.0/M_PI << setw(16) // radian
	 << endl;
    //  ofstream ofs9("dir-dat/fxx_fyy.dat");
    ofs9 << t   << setw(16)
	 << fxx << setw(16)
	 << fyy << endl;

    if (10*static_cast<int>(nt*0.1) == nt){
      FILE *gp5, *gp6, *gp7;
      sprintf(mfile,"dir-dat/s1-%05d.dat", nt); gp5 = fopen(mfile, "w");
      sprintf(mfile,"dir-dat/s2-%05d.dat", nt); gp6 = fopen(mfile, "w");
      sprintf(mfile,"dir-dat/s3-%05d.dat", nt); gp7 = fopen(mfile, "w");
      //ofstream ofs1("dir-dat/end_effector.dat");
      ofs1 << setw(16) << t
	   << setw(16) << fxx
	   << setw(16) << fyy
	   << setw(16) << 0.0
	   << setw(16) << r2.x1
	   << setw(16) << r2.y1
	   << setw(16) << 0.0
	   << setw(16) << keirosokudo
	   << endl;

      zu(gp5, gp6, gp7, env);
      fclose(gp5); fclose(gp6); fclose(gp7);
      zu_real(gp1, *(env.r1), *(env.r2), *(env.r3), *(env.th_));
      zu_virtual(gp2,*(env.v1),*(env.th_));
    }
  }
  fclose(gp1); fclose(gp2);

  printf("%s %d\n", __FUNCTION__, __LINE__);
  //  delete sc1;
}

double En_prism(const struct Return_sim_value_t &th_,
		Linear_curve_t& fk1, Linear_curve_t& fk2,
		Conf_value& ini, Conf_value& fin,
		Manipulator& v1,
		Manipulator& r1, Manipulator& r2)
{
  register double detJ, J11, J12, J21, J22, coef;
  register double cth1, cth12, cth2;
  register double sth1, sth12, sth2;
  double XXa, YYa, PL3;

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

  tau_1 = (a1+a2+a3+2.0*a4*cth2)*r1.tha + (a2+a3+a4*cth2)*r2.tha
    - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*sth2;
  tau_2 = (a2 + a3 + a4*cth2)*r1.tha
    + (a2 + a3)*r2.tha + a4*r1.thv*r1.thv*sth2;

  ev1 = b1*r1.thv + b2*r1.tha + b3*tau_1; ia1 = (ev1 - Kv*r1.thv)/Ra;
  ev2 = b1*r2.thv + b2*r2.tha + b3*tau_2; ia2 = (ev2 - Kv*r2.thv)/Ra;
  if (ev1*ia1 > 0.0){ Ene_1 = Ene_1 + ev1*ia1*dt; }
  if (ev2*ia2 > 0.0){ Ene_2 = Ene_2 + ev2*ia2*dt; }

  //Ene123 = Ene_1 + Ene_2 + Ene_3;
  Ene_123 = Ene_1 + Ene_2;

  return Ene_123;//Ene123;
}

double En_vertical(const struct Return_sim_value_t &z_,
		   BallScrew& bs, 
		   Manipulator &v1, Manipulator &r3)
{
  register double detJ, J11, J12, J21, J22, coef;
  register double cth1, cth12, cth2;
  register double sth1, sth12, sth2;  

  v1.z0 = z_.s;			// v1.z1, v1.z => zu_virtual
  
  //  r3.z0 = z_.s;
  r3.z1 = z_.s;
  r3.zv = z_.v;
  r3.za = z_.a;
  // r3.th, r3.thv, r3.tha is using motor value
  r3.th  = bs.frac*z_.s;
  r3.thv = bs.frac*z_.v;
  r3.tha = bs.frac*z_.a;

  //cout << bs.frac <<" "<<bs.pitch<< endl;
  //assert(false);

  
//  tau_1 = (a1+a2+a3+2.0*a4*cth2)*r1.tha + (a2+a3+a4*cth2)*r2.tha
//    - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*sth2;
//  tau_2 = (a2 + a3 + a4*cth2)*r1.tha
//    + (a2 + a3)*r2.tha + a4*r1.thv*r1.thv*sth2;
  tau_3 = bs.getTau(r3.tha);

  cout << bs.m <<"-"<< r3.m << endl;

  
  //  ev1 = b1*r1.thv + b2*r1.tha + b3*tau_1; ia1 = (ev1 - Kv*r1.thv)/Ra;
  //  ev2 = b1*r2.thv + b2*r2.tha + b3*tau_2; ia2 = (ev2 - Kv*r2.thv)/Ra;
  ev3 = b1*r3.thv + b2*r3.tha + b3*tau_3; ia3 = (ev3 - Kv*r3.thv)/Ra;
  //  if (ev1*ia1 > 0.0){ Ene_1 = Ene_1 + ev1*ia1*dt; }
  //  if (ev2*ia2 > 0.0){ Ene_2 = Ene_2 + ev2*ia2*dt; }
  if (ev3*ia3 > 0.0){ Ene_3 = Ene_3 + ev3*ia3*dt; }

  //Ene_123 = Ene_1 + Ene_2 + Ene_3;
  //  Ene_123 = Ene_1 + Ene_2 + Ene_3;
  Ene_123 = Ene_3;
  return Ene_123;//Ene_123;
}

//******************************************************************
void zu(FILE *_fp1, FILE *_fp2, FILE *_fp3, Environment &e)
{
  zu_virtual(_fp1, *(e.v1), *(e.th_));// virtual manipulator
  zu_real(_fp2, *(e.r1), *(e.r2), *(e.r3), *(e.th_));  // real manipulator
  zu_pathway(_fp3, *(e.v1));  // pathway
}

void zu_virtual(FILE *_fp, Manipulator& v1, Return_sim_value_t& th_)
{
 v1.th = th_.s;
 // v1.x1 = v1.x0 + v1.l*cos(v1.th);
 // v1.y1 = v1.y0 + v1.l*sin(v1.th);
 v1.z1 = v1.z = v1.z0;
 
 fprintf(_fp,"%lf %lf %lf\n", v1.x0, v1.y0, v1.z0);
 fprintf(_fp,"%lf %lf %lf\n", v1.x1, v1.y1, v1.z0);
 draw_circle_fp(_fp, v1.x0, v1.y0, v1.z0, RR);
 draw_circle_fp(_fp, v1.x1, v1.y1, v1.z0, RR*0.2); // end-effector
}

void zu_real(FILE *_fp,
	     Manipulator& r1,Manipulator& r2,Manipulator& r3,
	     Return_sim_value_t& th_)
{
 r1.x1 = r1.x0 + r1.l*cos(r1.th);
 r1.y1 = r1.y0 + r1.l*sin(r1.th);
 //  r1.z1 = r1.z = r1.z0;
 fprintf(_fp,"%lf %lf %lf\n", r1.x0, r1.y0, r1.z0);
 fprintf(_fp,"%lf %lf %lf\n", r1.x1, r1.y1, r1.z0);
 draw_circle_fp(_fp, r1.x0, r1.y0, r1.z0       , RR);
 draw_circle_fp(_fp, r1.x0, r1.y0, r1.z0-RR*0.1, RR);  
 
 r2.x0 = r1.x1; r2.x1 = r2.x0 + r2.l*cos(r1.th + r2.th);
 r2.y0 = r1.y1; r2.y1 = r2.y0 + r2.l*sin(r1.th + r2.th);
 //  r2.z1 = r2.z = r2.z0;
 fprintf(_fp," %lf %lf %lf\n", r2.x0, r2.y0, r2.z0);
 fprintf(_fp," %lf %lf %lf\n", r2.x1, r2.y1, r2.z1);  
 //  fprintf(fp2,"%lf %lf %lf\n", r2.x1, r2.y1, r2.z1);  
 draw_circle_fp(_fp, r2.x0, r2.y0, r2.z0      , RR*0.5);
 draw_circle_fp(_fp, r2.x0, r2.y0, r2.z0-RR*0.1, RR*0.5);  
 
 r3.x0 = r2.x1; r3.y0 = r2.y1; 
 r3.x1 = r3.x0; r3.y1 = r3.y0;
 r3.z0 = r3.z1 + r3.l*cos(0);
 r3.z  = r3.z1 + r3.r*cos(0);
 fprintf(_fp,"   %lf %lf %lf\n", r3.x0, r3.y0, r3.z0);
 fprintf(_fp,"   %lf %lf %lf\n", r3.x0, r3.y0, r3.z1);  
 //fprintf(_fp,"%lf %lf %lf\n", r3.x1, r3.y1, r3.z1);  
 draw_circle_fp(_fp, r3.x0, r3.y0, r3.z , RR*0.3);
 draw_circle_fp(_fp, r3.x0, r3.y0, r3.z0, RR*0.3);
 draw_circle_fp(_fp, r3.x0, r3.y0, r3.z1, RR*0.3);  
}
void zu_pathway(FILE *_fp, Manipulator& v1){
 draw_circle_fp(_fp, v1.x1, v1.y1, v1.z1, RR);
}

//-----------------------------------------------------------------
void draw_circle_fp(FILE *_fp, double x, double y, double z, double r)
{
  int i;
  double dth = M_PI/18.0, th = 0.0;
  fprintf(_fp,"\n");
  for(i = 0; i <= 37; i++){
    fprintf(_fp, "%lf %lf %lf\n", x + r*cos(th), y + r*sin(th), z);
    th = th + dth;
  }
  fprintf(_fp,"\n");  
}

void draw_circle_z
( double x, double y, double z, double r,
  double *XX, double *YY, double *ZZ)
{
  int i;
  double th = 0.0, dth = M_PI/180.0;
  for(i = 0; i <= 37; i++){
    XX[i] = x + r*cos(th); YY[i] = y + r*sin(th); ZZ[i] = z;
    th += dth;
  }
}


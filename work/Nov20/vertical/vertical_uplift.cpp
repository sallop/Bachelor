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

struct Linear_curve_t {double a, b;} fk1, fk2, fkv;
struct Conf_value{
  double x, y, z, vth, rth1, rth2;
} ini, fin;

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
} *sc1, *scv;
Sin_curve::Sin_curve(double _ttf,
		     double _th1_i,double _th1_f,
		     double _th1v_i, double _th1v_f)
  : ttf(_ttf), th1_i(_th1_i), th1_f(_th1_f),
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
	  , t, sin_th1(t), sin_th1v(t),sin_th1a(t) );
}

class BallScrew{
public:
  BallScrew(double _J, double _m, double _pitch);
  double J;			// motor side of inerita
  double m;			// mass of fuka side
  double frac;			// fraction of gear teeth  = 2.0*M_PI/pitch;
  double pitch;			// pitch of screw
  double getTau(double tha);		// motor side of torque
} *bs;
BallScrew::BallScrew(double _J, double _m, double _pitch)
  : J(_J), m(_m), pitch(_pitch)
{
  frac = 2.0*M_PI/pitch;
  std::cout << "J=" << J << setw(2) 
	    << "m=" << m << setw(2)
	    << "f=" << frac << endl;
}
double BallScrew::getTau(double tha)
{
  return (J + m/(frac*frac))*tha;
}

struct Manipulator{
  double l, r, m, Ig;
  double x0, y0, z0, x,  y,  z, x1, y1, z1;
  double xv1, yv1, xa1, ya1;
  double th, thv, tha;
  double lv, la;
  double zv, za;
  void print_value(FILE *_fp, char *name);
} v1, v2, r1, r2, r3;

void Manipulator::print_value(FILE *_fp, char *name)
{
  fprintf(_fp, "variable:%s\n", name);
  fprintf(_fp, "l=%lf,r=%lf,m=%lf,Ig=%lf\n", l, r, m, Ig);
  fprintf(_fp, "x0=%lf,y0=%lf,z0=%lf\n", x0, y0, z0);
  fprintf(_fp, "x  =%lf,y  =%lf,z  = %lf\n", x,  y,  z);
  fprintf(_fp, "x1 =%lf,y1 =%lf,z1 = %lf\n", x1, y1, z1);
  fprintf(_fp, "xv1=%lf,yv1=%lf\n",xv1, yv1);
  fprintf(_fp, "xa1=%lf,ya1=%lf\n", xa1, ya1);
  fprintf(_fp, "zv=%lf,za=%lf", zv, za);
  fprintf(_fp, "th=%lf,thv=%lf,tha=%lf\n", th, thv, tha);
  fprintf(_fp, "lv=%lf,la=%lf\n", lv, la);
}

double b1, b2, b3;		// moter
double a1, a2, a3, a4, a5;
double t, Ts, t1, ts[20], ttf, tt, Pdeg;
double th1_i, th1_f, th2_i, th2_f, fai_2i, fai_2f;
double th1v_i;
double ac1, Bc1, mL;
double omg1, Amp1, omg2, Amp2;	// using sincurve
double ev1, ev2, ev3, ia1, ia2, ia3;
// virtual manipulator horizontal
struct Return_sim_value_t {double s, v, a;} th_;
double Sth1[20], kiz1a, th_1, thv_1, tha_1, th1[20], th1v[20], th1a[20];
double Sth2[20], kiz2a, th_2, thv_2, tha_2, th2[20], th2v[20], th2a[20];
// virtual manipulator vertical
double Energ[15][15][15], Enemin2[50];
double Ene_1, Ene_2, Ene_3, enemin1, enemin2, ene_min2;
double Ene_123;

double tau_1, tau_2, tau_3, fz_3;// using En1, En2, En3
double cth1, cth2, cth12, sth1, sth2, sth12;
double detJ, XXa, YYa, PL3;

double fxx,fyy,fxa,fya,keirosokudo,keirosokudo2;

int n, nnn, i, nt, Mat[13][13][13], N3;
int sat[1+DIV], at[1+DIV];
int flag11, ff1, ff2, flag34;
int nt11, nt22, nt33;

void sinc();
double En1(const Return_sim_value_t &th_);
double En_vertical(const Return_sim_value_t &th_);

void zu(FILE *_fp1, FILE *_fp2, FILE *_fp3); void zu_real(FILE *_fp);
void zu_virtual(FILE *_fp); void zu_pathway(FILE *fp);
void draw_circle();
void draw_circle_fp(FILE *_fp, double x, double y, double z, double r);

//**********************************************************************
void Animation()
{
  FILE *gp;
  int i, k;
  gp = popen("gnuplot","w");
  
  fprintf(gp, "set view 45, 300\n");
  fprintf(gp,"set xrange[-0.15:0.15]\n");
  fprintf(gp,"set yrange[-0.07:0.20]\n");
  fprintf(gp,"set zrange[-0.15:0.15]\n");
  fprintf(gp,"set terminal x11\n");

  fprintf(gp,"set size 0.6, 0.6\n");
  fprintf(gp," set grid \n");
  for(k = 1;k <= 2; k++){
    for(i = 10; i <= 640; i += 10){    
      fprintf(gp,"splot ");
      fprintf(gp," 'dir-dat/s1-%05d.dat' w l 1",i);
      fprintf(gp,",'dir-dat/s2-%05d.dat' w l 3",i);
      fprintf(gp,",'dir-dat/s3-%05d.dat' w l 2",i);
      //      fprintf(gp,"0.10 5\n");
      fprintf(gp,"\n");
      fprintf(gp,"pause 0.05 \n");
    }
  }

  for(i = 10; i <= 640; i += 10){
    fprintf(gp, "splot ");
    fprintf(gp," 'dir-dat/s1-%05d.dat' w l 1",i);
    fprintf(gp,",'dir-dat/s2-%05d.dat' w l 3",i);
    fprintf(gp,",'dir-dat/s3-%05d.dat' w l 2",i);
    //    fprintf(gp,"0.10 5\n");
    fprintf(gp,"\n");    
    fprintf(gp,"set terminal gif\n");
    fprintf(gp,"set grid\n");
    fprintf(gp,"set output 'dir-gif/work-2_%05d.gif'\n",i);
    fprintf(gp,"rep\n");
  }
  
  for(k = 1;k <= 2; k++){
    for(i = 10; i <= 640; i += 10){      
      fprintf(gp,"splot ");
      fprintf(gp," 'dir-dat/s1-%05d.dat' w l 1",i);
      fprintf(gp,",'dir-dat/s2-%05d.dat' w l 3",i);
      fprintf(gp,",'dir-dat/s3-%05d.dat' w l 2",i);
      //      fprintf(gp,"0.10 5\n");
      fprintf(gp,"\n");      
      fprintf(gp,"pause 0.05 \n");
    }
  }

  for(i = 10; i <= 640; i += 10){
    fprintf(gp, "splot ");
    fprintf(gp, " 'dir-dat/s1-%05d.dat' w l 1",i);
    fprintf(gp, ",'dir-dat/s2-%05d.dat' w l 3",i);
    fprintf(gp, ",'dir-dat/s3-%05d.dat' w l 2",i);
    //    fprintf(gp, "0.10 5\n");
    fprintf(gp, "\n");    
    fprintf(gp,"set terminal gif\n");
    fprintf(gp,"set grid\n");
    fprintf(gp,"set output 'dir-gif/work-2.gif'\n");
    fprintf(gp,"rep\n");
  }

  pclose(gp);  
}

inline double
det(double a11, double a12, double a21, double a22){
  return a11*a22 - a12*a21;
}

void init_variable(FILE *_fp)
{
  //-initialize value
  mL = 0.0; //g = 9.8;
  // setting time variable
  //Ts = 0.64 ; t1=(Ts/(double)DIV); //t1 = Ts/16.0;
  Ts = 0.64 ;
  t1 = Ts/static_cast<double>(DIV); //t1 = Ts/16.0;
  fprintf(_fp,"%d Ts=%lf,t1=%lf,DIV=%d\n",__LINE__, Ts, t1, DIV);
  for(i = 1; i <= DIV; i++){
    ts[i] = t1*i;
    fprintf(_fp,"%d ts[%d]=%lf\n", __LINE__, i, ts[i]);
  }
  ttf = ts[DIV];
  
  // actuator
  //b1 = Kv + Ra*Dm/Kt; b3 = Ra/Kt; b2 = b3*Im;// dt=0.001;
  b1 = Kv + Ra*Dm/Kt; b2 = Ra*Im/Kt; b3 = Ra/Kt; // dt=0.001;
  fprintf(_fp,"%d b1=%lf,b2=%lf,b3=%lf\n",__LINE__,b1,b2,b3);
  
  v1.x0 = -0.05 - 0.005; v1.y0 = 0.10 + 0.05; v1.z0 = 0.0;    
  //  v1.x0 = -0.05 - 0.01; v1.y0 = 0.10 + 0.05; v1.z0 = 0.0;
  //  v1.x0 = -0.05; v1.y0 = 0.10 + 0.05; v1.z0 = 0.0;

  fprintf(_fp, "%d v1.x0=%lf,v1.y0=%lf,v1.z0=%lf\n"
	  ,__LINE__, v1.x0, v1.y0, v1.z0);
  r1.m = m1; r2.m = m2; r3.m = m3;
  r1.l = L1; r1.r = r1.l/2.0; r1.Ig = r1.m*r1.l*r1.l/12.0;
  r2.l = L2; r2.r = r2.l/2.0; r2.Ig = r2.m*r2.l*r2.l/12.0;
  r3.l = L3; r3.r = r3.l/2.0; r3.Ig = 0.0;
  // start and end point
  {
    double l_hypot, l_base;
    ini.x  =  0.05; ini.y  = 0.10; ini.z  =  0.00;
    // using work.cpp (fxx, fyy)
    fin.x  = -0.05; fin.y  = 0.10; fin.z  =  0.0618;

    //l_hypot = sqrt(SQR(ini.xt - v1.x0) + SQR(ini.yt  - v1.y0));
    l_hypot = hypot(ini.x - v1.x0, ini.y  - v1.y0);
    //l_base  = sqrt(SQR(fin.xal - v1.x0)+ SQR(fin.yal - v1.y0));
    l_base  = hypot(fin.x - v1.x0, fin.y - v1.y0);
    fprintf(_fp,"%d hypot=%lf,base=%lf,l_base/l_hypot=%lf\n"
	    ,__LINE__, l_hypot, l_base, l_base/l_hypot);

    //    cfg_vth.init = acos(l_base/l_hypot);
    ini.vth = atan2(ini.y - v1.y0, ini.x - v1.x0);
    fin.vth = atan2(fin.y - v1.y0, fin.x - v1.x0);

    fprintf(_fp,"%d ini.vtht=%lf,(d)%lf\n"
	    ,__LINE__, ini.vth, ini.vth*180.0/M_PI);
    fprintf(stdout,"%d ini.vtht=%lf,(d)%lf\n"
	    ,__LINE__, ini.vth, ini.vth*180.0/M_PI);
    fprintf(_fp,"%d fin.vth=%lf,(d)%lf\n"
	    ,__LINE__, fin.vth, fin.vth*180.0/M_PI);
    fprintf(stdout,"%d fin.vth=%lf,(d)%lf\n"
	    ,__LINE__, fin.vth, fin.vth*180.0/M_PI);    
    l_hypot = hypot(ini.x, ini.y);
    fprintf(_fp,"%d hypot(%lf,%lf)=%lf,/2=%lf\n"
	    ,__LINE__, ini.x, ini.y, l_hypot, l_hypot/2.0);
    ini.rth1 = atan2(ini.y, ini.x) - acos(0.5*l_hypot/r1.l);
    fprintf(_fp,"%d rth1.init=%lf,(d)%lf\n"
	    ,__LINE__, ini.rth1, ini.rth1*180.0/M_PI);
    ini.rth2 = 2.0*acos(0.5*l_hypot/r1.l);
    fprintf(_fp,"%d rth2.init=%lf,(d)%lf\n"
	    ,__LINE__, ini.rth2, ini.rth2*180.0/M_PI);
  }

  th_.s = ini.vth;
  // virtual manipulator
  v1.m = m1; v1.th = th_.s;
  fprintf(_fp,"%d v1.m=%lf,v1.th=%lf\n",__LINE__, v1.m, v1.th);
  v1.l = 0.05; v1.r = v1.l/2.0;
  fprintf(_fp, "%d v1.l=%lf,v1.r=%lf\n",__LINE__, v1.l, v1.r);

  v2 = v1;

  r1.x0 = 0.0;
  r1.x  = r1.x0 + r1.r*cos(ini.rth1);
  r1.x1 = r1.x0 + r1.l*cos(ini.rth1);
  r1.y0 = 0.0;
  r1.y  = r1.y0 + r1.r*sin(ini.rth1);
  r1.y1 = r1.y0 + r1.l*sin(ini.rth1);
  r1.z0 = HH; r1.z = HH; r1.z1 = HH;

  r2.x0 = r1.x1;
  r2.x  = r2.x0 + r2.r*cos(ini.rth1 + ini.rth2);
  r2.x1 = r2.x0 + r2.l*cos(ini.rth1 + ini.rth2);
  r2.y0 = r1.y1;
  r2.y  = r2.y0 + r2.r*sin(ini.rth1 + ini.rth2);
  r2.y1 = r2.y0 + r2.l*sin(ini.rth1 + ini.rth2);  
  r2.z0 = HH; r2.z = HH; r2.z1 = HH;

  r3.x0 = r2.x1; r3.x  = r2.x1; r3.x1 = r2.x1;
  r3.y0 = r2.y1; r3.y  = r2.y1; r3.y1 = r2.y1;
  r3.z0 = HH   ; r3.z  = r3.z0 - r3.l/2.0; r3.z1 = r3.z0 - r3.l;

  fk1.a = 0.0; fk1.b = 0.1;
  fk2.a = tan(th_.s); fk2.b = v1.y0 - fk2.a*v1.x0;
  fprintf(_fp,"%d fk1.a=%lf,fk1.b=%lf\n", __LINE__, fk1.a, fk1.b);
  fprintf(_fp,"%d fk2.a=%lf,fk2.b=%lf\n", __LINE__, fk2.a, fk2.b);
  fprintf(_fp,"%d th_.s=%8.4f tan(th_.s)=%8.4f\n"
	  , __LINE__, th_.s, tan(th_.s));

  fprintf(_fp,"%d cfg_rthi=%lf,fth2_i=%lf\n"
	  , __LINE__, ini.rth1, ini.rth2);  
  
  N3 = (int)(NN*0.5);
  fprintf(_fp,"%d NN=%3d  N3=%3d\n",__LINE__,NN,N3);

  v1.print_value(_fp,"v1"); r1.print_value(_fp,"r1");
  r2.print_value(_fp,"r2"); r3.print_value(_fp,"r3");
  
  //------------- Scara Robot --------------------------
  a1 = r1.Ig + r1.m*r1.r*r1.r + r2.m*r1.l*r1.l + r3.m*r1.l*r1.l;
  a2 = r2.Ig + r2.m*r2.r*r2.r + r3.m*r2.l*r2.l;
  a3 = r3.Ig;
  a4 = r2.m*r1.l*r2.r + r3.m*r1.l*r2.l;
  //a5 = ;
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
  //-----------------------------------------------------
  th1_i = ini.vth; th1_f = fin.vth;
  fprintf(_fp,"%d th1_i=%lf, th1_f=%lf\n",__LINE__, th1_i, th1_f);
  fprintf(_fp,"%d th1_i(d)=%lf, th1_f(d)=%lf\n",__LINE__
	  , th1_i*180.0/M_PI, th1_f*180.0/M_PI);

  th1v_i = 0.0 ;
  fprintf(_fp,"%d Amp1=%lf, th1v_i=%lf\n",__LINE__, Amp1, th1v_i);
  th1[0]    = th1_i ; th1v[0]   = th1v_i;
  th1[DIV]  = th1_f ; th1v[DIV] = 0.0;

  //  nt11 = 214; nt22 = 427;// Sep/5
  nt11 = 220; fprintf(_fp,"%d nt22=%d\n",__LINE__, nt22);
  nt22 = 440; fprintf(_fp,"%d nt11=%d\n",__LINE__, nt11);
  
  // libration and amplitude
  //max = (int)(Ts/dt) ;
  Pdeg = 180.0/M_PI ;
  omg1 = M_PI/ttf;
  //fprintf(_fp,"%d max=%d,Pdeg=%lf,omg1=%lf\n",__LINE__,max, Pdeg, omg1);
  fprintf(_fp,"%d Pdeg=%lf,omg1=%lf\n",__LINE__, Pdeg, omg1);
  Amp1 = ((th1_f - th1_i)/2.0)*omg1;
  //Amp1 = ((fin.vth - ini.vth)/2.0)*omg1;
  ac1 = (Amp1 - th1v_i)/ttf;
  Bc1 = (M_PI/2.0)*((th1_f - th1_i)/ttf - (Amp1 - th1v_i)/2.0);
  fprintf(_fp,"%d ac1=%lf,ttf=%lf,Bc1=%lf\n",__LINE__,ac1,ttf,Bc1);

  sc1 = new Sin_curve( ttf, ini.vth, fin.vth, 0.0, 0.0); 
  // ttf(:time to final), (th1), (th1v)
  scv = new Sin_curve( ttf, ini.z  , fin.z  , 0.0, 0.0);
  bs = new BallScrew(Im, r1.m, 2.0*M_PI);
 // motor side: moment of inerita 
  
  for(n = 1;n <= DIV; n++){
    Sth1[n] = sc1->sin_th1(ts[n]);
    sc1->print_value(stdout, ts[n]);
    fprintf(_fp,"%d tt=%lf,Sth1[%d]=%lf %lf\n"
	    ,__LINE__, tt, n, Sth1[n], ts[n]);
  }
  kiz1a = 0.04; fprintf(_fp,"%d kiz1a=%lf\n",__LINE__,kiz1a);
  kiz2a = 0.04; fprintf(_fp,"%d kiz2a=%lf\n",__LINE__,kiz2a);
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
  FILE *fp1;
  fp1 = fopen("work2_variable.dat","w");

  start = gettimeofday_sec();
  init_variable(fp1);
  fclose(fp1);
  sinc();
  end = gettimeofday_sec();
  Animation();
  
  printf("time=%g,start=>%g,end=%g\n", end - start, start, end);
  return 0 ;
}

//************************************************************
void test_sinc(FILE *_fp)
{
  fprintf(_fp,"sin_th1=%lf\n", th_.s);
  fprintf(_fp,"sin_th1v=%lf\n", th_.v);
  fprintf(_fp,"sin_th1a=%lf\n", th_.a);
  fprintf(_fp,"Amp1=%lf\n", Amp1);
  fprintf(_fp,"t=%lf\n", t);
  fprintf(_fp,"ttf=%lf\n", ttf);
  fprintf(_fp,"th1_i=%f\n", th1_i);
  fprintf(_fp,"th1v_i=%lf\n", th1v_i);
  fprintf(_fp,"th1_f=%f\n", th1_f);  
}

void sinc()
{
  char mfile[32];
  FILE *gp1, *gp2;
  
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
  r1.th  = ini.rth1;  r1.thv = 0.0;
  r2.th  = ini.rth2;  r2.thv = 0.0;

  //  test_sinc(stdout);
  for(nt = 0; nt <= static_cast<int>(ttf/dt); nt++){
    t = t + dt;       
//    th_.s = sc1->sin_th1(t);
//    th_.v = sc1->sin_th1v(t);
//    th_.a = sc1->sin_th1a(t);
    th_.s = scv->sin_th1(t);
    th_.v = scv->sin_th1v(t);
    th_.a = scv->sin_th1a(t);
    //    En1(th_);
    En_vertical(th_);
    
    ofs0 << setw(12) << t      
	 << setw(12) << r1.th  << setw(12) << r2.th 
	 << setw(12) << r1.thv << setw(12) << r2.thv
	 << setw(12) << r1.tha << setw(12) << r2.tha
	 << endl;

    ofs_ene1 << t   << setw(16)  << Ene_1 << setw(16)
	     << ev1 << setw(16)  << ia1   << setw(16)
	     << tau_1 << Ene_123 << endl;

    ofs_ene2 << t   << setw(16)  << Ene_2 << setw(16)
	     << ev2 << setw(16)  << ia2   << setw(16)
	     << tau_2 << Ene_123 << endl;

    ofs_ene3 << t   << setw(16)  << Ene_3 << setw(16)
	     << ev3 << setw(16)  << ia3   << setw(16)
	     << tau_3 << Ene_123 << endl;

    ofs4 << t << setw(16)
	 << tan(th_.s) << setw(16)
	 << th_.s*180.0/M_PI << setw(16)
	 << th_.v << setw(16)
	 << th_.a << endl;

    ofs8 << t << setw(16)
	 << th_.s << setw(16)
	 << th_.v << setw(16)
	 << th_.a << setw(16)
	 << th_.s*180.0/M_PI << setw(16) // radian
	 << endl;

    ofs9 << t   << setw(16)
	 << fxx << setw(16)
	 << fyy << endl;

    if (10*(int)(nt*0.1) == nt){
      FILE *gp5, *gp6, *gp7;
      sprintf(mfile,"dir-dat/s1-%05d.dat", nt); gp5 = fopen(mfile, "w");
      sprintf(mfile,"dir-dat/s2-%05d.dat", nt); gp6 = fopen(mfile, "w");
      sprintf(mfile,"dir-dat/s3-%05d.dat", nt); gp7 = fopen(mfile, "w");
      ofs1 << setw(16) << t
	   << setw(16) << fxx
	   << setw(16) << fyy
	   << setw(16) << 0.0
	   << setw(16) << r2.x1
	   << setw(16) << r2.y1
	   << setw(16) << 0.0
	   << setw(16) << keirosokudo
	   << endl;

      zu(gp5, gp6, gp7);
      fclose(gp5); fclose(gp6); fclose(gp7);
      zu_real(gp1); zu_virtual(gp2);
    }
  }
  fclose(gp1); fclose(gp2);

  printf("%s %d\n", __FUNCTION__, __LINE__);
  //  delete sc1;
}

// **************************************************
double En1(const struct Return_sim_value_t &th_)
{
  register double detJ, J11, J12, J21, J22, coef;
  register double cth1, cth12, cth2;
  register double sth1, sth12, sth2;  


  cth1 = cos(r1.th); cth12 = cos(r1.th + r2.th); cth2 = cos(r2.th);
  sth1 = sin(r1.th); sth12 = sin(r1.th + r2.th); sth2 = sin(r2.th);

  fk1.a = (fin.y - ini.y)/(fin.x - ini.x);
  fk1.b = ini.y;
  fk2.a = tan(th_.s);
  fk2.b = v1.y0 - fk2.a*v1.x0 ;
  fxx   = (fk2.b - fk1.b)/(fk1.a - fk2.a);
  fyy   = fk2.a*fxx + fk2.b;

  v1.l   = hypot(fxx-v1.x0 , fyy-v1.y0);
  coef   = (fk1.a*tan(th_.s) + 1.0)/(fk1.a - tan(th_.s));
  v1.lv  = coef*v1.l*th_.v;
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

  //Ene_123 = Ene_1 + Ene_2 + Ene_3;
  Ene_123 = Ene_1 + Ene_2;

  return Ene_123;//Ene_123;
}

double En_prism(const struct Return_sim_value_t &th_)
{
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

  tau_1 = (a1+a2+a3+2.0*a4*cth2)*r1.tha + (a2+a3+a4*cth2)*r2.tha
    - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*sth2;
  tau_2 = (a2 + a3 + a4*cth2)*r1.tha
    + (a2 + a3)*r2.tha + a4*r1.thv*r1.thv*sth2;

  ev1 = b1*r1.thv + b2*r1.tha + b3*tau_1; ia1 = (ev1 - Kv*r1.thv)/Ra;
  ev2 = b1*r2.thv + b2*r2.tha + b3*tau_2; ia2 = (ev2 - Kv*r2.thv)/Ra;
  if (ev1*ia1 > 0.0){ Ene_1 = Ene_1 + ev1*ia1*dt; }
  if (ev2*ia2 > 0.0){ Ene_2 = Ene_2 + ev2*ia2*dt; }

  //Ene123 = Ene_1 + Ene_2 + Ene_3;
  Ene123 = Ene_1 + Ene_2;

  return Ene123;//Ene123;
}

double En_vertical(const struct Return_sim_value_t &z_)
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
  r3.th  = bs->frac*z_.s;
  r3.thv = bs->frac*z_.v;
  r3.tha = bs->frac*z_.a;
  
  tau_1 = (a1+a2+a3+2.0*a4*cth2)*r1.tha + (a2+a3+a4*cth2)*r2.tha
    - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*sth2;
  tau_2 = (a2 + a3 + a4*cth2)*r1.tha
    + (a2 + a3)*r2.tha + a4*r1.thv*r1.thv*sth2;
  tau_3 = bs->getTau( bs->frac*z_.a );
  
  ev1 = b1*r1.thv + b2*r1.tha + b3*tau_1; ia1 = (ev1 - Kv*r1.thv)/Ra;
  ev2 = b1*r2.thv + b2*r2.tha + b3*tau_2; ia2 = (ev2 - Kv*r2.thv)/Ra;
  ev3 = b1*r3.thv + b2*r3.tha + b3*tau_3; ia3 = (ev3 - Kv*r3.thv)/Ra;
  if (ev1*ia1 > 0.0){ Ene_1 = Ene_1 + ev1*ia1*dt; }
  if (ev2*ia2 > 0.0){ Ene_2 = Ene_2 + ev2*ia2*dt; }
  if (ev3*ia3 > 0.0){ Ene_3 = Ene_3 + ev3*ia3*dt; }

  //Ene_123 = Ene_1 + Ene_2 + Ene_3;
  Ene_123 = Ene_1 + Ene_2 + Ene_3;

  return Ene_123;//Ene_123;
}

//******************************************************************
void zu(FILE *_fp1, FILE *_fp2, FILE *_fp3){
  zu_virtual(_fp1);// virtual manipulator
  zu_real(_fp2);  // real manipulator
  zu_pathway(_fp3);  // pathway
}
void zu_virtual(FILE *_fp){
  v1.th = th_.s;
  v1.x1 = v1.x0 + v1.l*cos(v1.th);
  v1.y1 = v1.y0 + v1.l*sin(v1.th);
  v1.z1 = v1.z = v1.z0;
  
  fprintf(_fp,"%lf %lf %lf\n", v1.x0, v1.y0, v1.z0);
  fprintf(_fp,"%lf %lf %lf\n", v1.x1, v1.y1, v1.z0);
  //  fprintf(_fp,"%lf %lf %lf\n", v1.x1, v1.y1, v1.z1);
  draw_circle_fp(_fp, v1.x0, v1.y0, v1.z0, RR);
  draw_circle_fp(_fp, v1.x1, v1.y1, v1.z0, RR*0.2);  
}
void zu_real(FILE *_fp){
  r1.x1 = r1.x0 + r1.l*cos(r1.th);
  r1.y1 = r1.y0 + r1.l*sin(r1.th);
  //  r1.z1 = r1.z = r1.z0;
  fprintf(_fp,"%lf %lf %lf\n", r1.x0, r1.y0, r1.z0);
  fprintf(_fp,"%lf %lf %lf\n", r1.x1, r1.y1, r1.z0);
  //fprintf(_fp,"%lf %lf %lf\n", r1.x1, r1.y1, r1.z1);
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
void zu_pathway(FILE *_fp){
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

#include <cstdio>
#include <cstring>
#include <string>
#include <cfloat>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <fenv.h>
#include <unistd.h>
#include "GraphicalPrimitive.hpp"
#include "Manipulator.hpp"
#include "Energy.hpp"
#include "util.h"
#include "variables.hpp"
#include "SinCurv.hpp"

using namespace std;

void make_gif_file();
int plot_flags[1+60];

enum PlotFlagType{ ePLOT=0, eSPLOT, eGif, eSize };
struct PlotFlag
{
  int flag[eSize];
  PlotFlag()
  {
    for(uint i=0; i < eSize; ++i)
      {
	flag[i] = 0;
      }
  }
  void set_flag( enum PlotFlagType pflag)
  {
    flag[pflag] = 1;
  }
} plotflag;
int viewflag = 0;

double th1_i, th1_f, th1v_i, th1v_f, th1a_i, th1a_f;
double th2_i, th2_f, th2v_i, th2v_f, th2a_i, th2a_f;
double th3_i, th3_f, th3v_i, th3v_f, th3a_i, th3a_f;

//Model model;
  
void view_animation( std::string prefix )
{
  struct {
    double x, y;
  } static const view[] = {
    {  67.0, 133.0 },
    {  90.0, 360.0 },
    { 107.0, 277.0 },
    {   0.0,   0.0 },
    {   0.0,  90.0 },
  };

  FILE *gp = popen("gnuplot","w");
  fprintf(gp, "set xrange[%lf:%lf]\n", -0.08, 0.10);
  fprintf(gp, "set yrange[%lf:%lf]\n", -0.02, 0.16);
  fprintf(gp, "set zrange[%lf:%lf]\n",  0.00, 0.20);
  //fprintf(gp, "set view %lf, %lf\n"  , view[0].x, view[0].y);
  //fprintf(gp, "set view %lf, %lf\n"  , view[1].x, view[1].y);
  //fprintf(gp, "set view %lf, %lf\n"  , view[2].x, view[2].y);
  //fprintf(gp, "set view %lf, %lf\n"  , view[3].x, view[3].y);
  fprintf(gp, "set view %lf, %lf\n"  , view[viewflag].x, view[viewflag].y);
  fprintf(gp, "set grid\n");

  //  fprintf(gp, "set terminal gif\n");

  if( plotflag.flag[ePLOT] )
    {
      for(int ii=10; ii < static_cast<int>(Ts/dt); ii+=10)
	{ // Magick
	  fprintf(gp, "plot ");
	  fprintf(gp, "'dir-dat/%s_sinc_e1_%05d.dat' w l", prefix.c_str(), ii);
	  fprintf(gp, ",");    
	  fprintf(gp, "'dir-dat/%s_sinc_e2_%05d.dat' w l", prefix.c_str(), ii);
	  fprintf(gp, ",");    
	  fprintf(gp, "'dir-dat/%s_sinc_e3_%05d.dat'", prefix.c_str(), ii);
	  fprintf(gp, "\n");    
	  //    fprintf(gp, "pause 0.01\n");
	  cout << __FUNCTION__<<"-"<< prefix <<"-"<< ii << endl;
	}

      for(int ii=10; ii < static_cast<int>(Ts/dt); ii+=10)
	{
	  //for(int ii=10; ii < static_cast<int>(3.0*ts[8]/dt); ii+=10){//Magick
	  //fprintf(gp, "set output 'dir-gif/%s-%05d.gif\n", prefix.c_str(), ii);
	  fprintf(gp, "plot ");
	  fprintf(gp, "'dir-dat/%s_saiteki_e1_%05d.dat' w l", prefix.c_str(), ii);
	  fprintf(gp, ",");
	  fprintf(gp, "'dir-dat/%s_saiteki_e2_%05d.dat' w l ls 3", prefix.c_str(), ii);
	  fprintf(gp, ",");    
	  fprintf(gp, "'dir-dat/%s_saiteki_e3_%05d.dat'", prefix.c_str(), ii);
	  fprintf(gp, "\n");
	  //    fprintf(gp, "pause 0.01\n");
	  cout << __FUNCTION__<<"-"<< prefix <<"-"<< ii << endl;
	}
    }

  if ( plotflag.flag[eSPLOT] )
    {
      for(int ii=10; ii < static_cast<int>(Ts/dt); ii+=10)
	{ // Magick
	  fprintf(gp, "splot ");
	  fprintf(gp, "'dir-dat/%s_sinc_e1_%05d.dat' w l", prefix.c_str(), ii);
	  fprintf(gp, ",");    
	  fprintf(gp, "'dir-dat/%s_sinc_e2_%05d.dat' w l", prefix.c_str(), ii);
	  fprintf(gp, ",");    
	  fprintf(gp, "'dir-dat/%s_sinc_e3_%05d.dat'", prefix.c_str(), ii);
	  fprintf(gp, "\n");    
	  //    fprintf(gp, "pause 0.01\n");
	  cout << __FUNCTION__<<"-"<< prefix <<"-"<< ii << endl;
	}
      
      for(int ii=10; ii < static_cast<int>(Ts/dt); ii+=10)
	{
	  //for(int ii=10; ii < static_cast<int>(3.0*ts[8]/dt); ii+=10){//Magick
	  //fprintf(gp, "set output 'dir-gif/%s-%05d.gif\n", prefix.c_str(), ii);
	  fprintf(gp, "splot ");
	  fprintf(gp, "'dir-dat/%s_saiteki_e1_%05d.dat' w l", prefix.c_str(), ii);
	  fprintf(gp, ",");
	  fprintf(gp, "'dir-dat/%s_saiteki_e2_%05d.dat' w l ls 3", prefix.c_str(), ii);
	  fprintf(gp, ",");    
	  fprintf(gp, "'dir-dat/%s_saiteki_e3_%05d.dat'", prefix.c_str(), ii);
	  fprintf(gp, "\n");
	  //    fprintf(gp, "pause 0.01\n");
	  cout << __FUNCTION__<<"-"<< prefix <<"-"<< ii << endl;
	}
    }

  if( plotflag.flag[eGif] )
    {
      for(int ii=10; ii < static_cast<int>(Ts/dt); ii+=10)
	{ // Magick
	  //fprintf(gp, "splot ");
	  fprintf(gp, "set terminal gif\n");
	  fprintf(gp, "set output 'dir-gif/%s_sinc_eN_%05d.dat'\n", prefix.c_str(), ii);
	  fprintf(gp, "plot ");
	  fprintf(gp, "'dir-dat/%s_sinc_e1_%05d.dat' w l", prefix.c_str(), ii);
	  fprintf(gp, ",");    
	  fprintf(gp, "'dir-dat/%s_sinc_e2_%05d.dat' w l", prefix.c_str(), ii);
	  fprintf(gp, ",");    
	  fprintf(gp, "'dir-dat/%s_sinc_e3_%05d.dat' w l", prefix.c_str(), ii);
	  fprintf(gp, "\n");    
	  //    fprintf(gp, "pause 0.01\n");
	  cout << __FUNCTION__<<"-"<< prefix <<"-"<< ii << endl;
	}
    }
  

  

  pclose(gp);
  cout << __FUNCTION__<<"-"<< prefix << endl;
}

void init()
{
  double tha_1, thv_1, th_1;	// return sim() parametar
  double tha_2, thv_2, th_2;	// return sim() parametar
  double tha_3, thv_3, th_3;	// return sim() parametar
  
  //t  = 0.0;
  //  Ts = 0.12; t1 = Ts/8.0;// pull pulley_test2.cpp
  Ts = 0.64; t1 = Ts/8.0;
  for(int i=0; i <= 8; i++){ ts[i] = t1*i; }

  kiza1 = 0.01; kiza2 = 0.01; kiza3 = 0.01;
  // Revolute(mass, length, x0, y0, z0, rth);
  v1 = VRManipulator( 0.0, line_i_v0,
		      pnt_v0.x, pnt_v0.y, pnt_v0.z,
		      vth_i );
  v1.setXYZ(pnt_i.x, pnt_i.y, pnt_i.z);

  //  v1.setXYZ(pnt_i.x, pnt_i.y, pnt_i.z);
  // pnt_r0 is a joint base position
  r1 = Revolute(m1, L1, pnt_r0.x, pnt_r0.y, pnt_r0.z, rth1);
  r1.setXYZ(pnt_r1.x, pnt_r1.y, pnt_r1.z);
  r2 = Revolute(m2, L2, pnt_r1.x, pnt_r1.y, pnt_r1.z, rth2);
  r2.setXYZ(pnt_i.x, pnt_i.y, pnt_i.z);
  r3 = Prismatic(m3, L3, pnt_i.x, pnt_i.y, pnt_i.z, radius);
  r3.setXYZ(pnt_i.x, pnt_i.y, pnt_i.z + r3.l);
  
  sc1 = SinCurv( Ts, cfg_th1.ini, cfg_th1.fin, 0.0, 0.0 );// up
  sc2 = SinCurv( Ts, cfg_th2.ini, cfg_th2.fin, 0.0, 0.0 ); // prismatic
  sc3 = SinCurv( Ts, cfg_th3.ini, cfg_th3.fin, 0.0, 0.0 );// down

  th_1  = cfg_th1.ini      ; // initialize value
  thv_1 = sc1.sin_th1v(0.0); // 
  tha_1 = sc1.sin_th1a(0.0); // 
  th_2  = cfg_th2.ini      ; // initialize value
  thv_2 = sc2.sin_th1v(0.0); // 
  tha_2 = sc2.sin_th1a(0.0); // 
  th_3  = cfg_th3.ini      ; // initialize value
  thv_3 = sc3.sin_th1v(0.0); // 
  tha_3 = sc3.sin_th1a(0.0); // 
  
  {// initialize theta
    //    th1_i  = vth_i; th1_f  = vth_f ;
    th1_i = cfg_th1.ini; th1_f = cfg_th1.fin;
    th1v_i = 0.0       ; th1v_f = 0.0       ;
    th1[0] = th1_i     ; th1[8] = th1_f     ;// th1[DIV] = th1_f;
    th1v[0] = 0.0      ; th1v[8] = th1v_f   ; //    th1v[DIV] = th1v_f;
    // th2[]
    th2_i = cfg_th2.ini; th2_f = cfg_th2.fin;
    th2v_i = 0.0       ; th2v_f = 0.0   ;
    th2[0] = th2_i     ; th2[8] = th2_f;// th2[DIV] = th2_f;
    th2v[0] = 0.0      ; th2v[8] = th2v_f; //    th2v[DIV] = th2v_f;
    // th3[]
    th3_i = cfg_th3.ini; th3_f = cfg_th3.fin;
    th3v_i = 0.0       ; th3v_f = 0.0   ;
    th3[0] = th3_i     ; th3[8] = th3_f; // th3[DIV] = th3_f;
    th3v[0] = 0.0      ; th3v[8] = th3v_f;  // th3v[DIV] = th3v_f;
  }

  //fk1.a = (pnt_f.y - pnt_i.y)/(pnt_f.x - pnt_i.x);// direct
  //fk1.b = pnt_i.y;				  // direct
  fk1.a = (pnt_c.y - pnt_i.y)/(pnt_c.x - pnt_i.x); // inclination
  fk1.b =  pnt_i.y;		// intercept
  fk2.a = tan(v1.th);		// fk2.a = tan(th_.s); variable
  fk2.b = v1.y0 - fk2.a*v1.x0 ; // => check x0         variable
  fk3.a = (pnt_f.y - pnt_c.y)/(pnt_f.x - pnt_c.x);
  fk3.b =  pnt_c.y;
  fk4.a = (pnt_f.y - pnt_i.y)/(pnt_f.x - pnt_i.x);// direct
  fk4.b =  pnt_i.y;				  // direct
  
  {//init coefficient
    a1 = r1.Ig + r1.m*r1.r*r1.r + r2.m*r1.l*r1.l + r3.m*r1.l*r1.l;
    a2 = r2.Ig + r2.m*r2.r*r2.r + r3.m*r2.l*r2.l;
    a3 = r3.Iz;
    a4 = r2.m*r1.l*r2.r + r3.m*r1.l*r2.l;   
  }
  
  ofstream ofs1("init_position.dat");
  //  zu(ofs1, set_position_pass);

  //model = Model(v1, r1, r2, r3, e1, e2, e3);

  
  //  if ( Ts1/dt < 1000 || Ts2/dt < 1000 || Ts3/dt < 1000){    
  //if ( (Ts/dt) < 1000 ){
  // 
  //  printf("sample point not enough\n");
  //  assert(false);
  //}

  ofstream ofs("init_end.dat");
  // zu(ofs, set_position_pass );
  //  assert(false);

  for(uint i = 0; i < sizeof(plot_flags)/sizeof(*plot_flags); ++i)
    plot_flags[i] = 0;

  plot_flags[ 1] = 1;
  plot_flags[10] = 1;
  plot_flags[20] = 1;
  plot_flags[30] = 1;
  plot_flags[40] = 1;
  plot_flags[50] = 1;
  plot_flags[60] = 1;  
}

int
main(int argc, char **argv)
{
  feenableexcept (FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  
  struct { double start, end; } bench1, bench2;
  Manipulator m(m1, LL, 0.0, 0.0, 0.0);
  int opt;
  while ( (opt = getopt(argc, argv, "psgv:")) != -1)
    {
      switch(opt)
	{
	case 'p':
	  plotflag.set_flag( ePLOT );
	  break;
	case 's':
	  plotflag.set_flag( eSPLOT );
	case 'g':
	  plotflag.set_flag( eGif );
	  break;
	case 'v':
	  viewflag = atoi(optarg);
	  break;
	default:
	  break;
	}
    }

  argv += optind;
  
  init();
  
  view_animation(argv[0]);
  //view_animation("nagato");
  
  return 0;
}

void make_gif_file()
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
	if (plot_flags[jj])
	  fprintf(gp, ",'dir-dat/saiteki_%05d.dat' u 1:%d w l lt 2"
		  , jj, data[ii].col);

      fprintf(gp, ",'dir-dat/saiteki_%05d.dat' u 1:%d w l lw 2 lt 3"
	      , REPEAT, data[ii].col);
      fprintf(gp, ",0  w l lw 1 lt 55");      
      fprintf(gp, "\n");
    }
  
  pclose(gp);
}

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
#include "variables.hpp"

using namespace std;

void make_gif_file();
int plot_flags[1+60];

double th1_i, th1_f, th1v_i, th1v_f, th1a_i, th1a_f;
double th2_i, th2_f, th2v_i, th2v_f, th2a_i, th2a_f;
double th3_i, th3_f, th3v_i, th3v_f, th3a_i, th3a_f;

  
void view_animation( std::string prefix )
{
  struct {
    double x, y;
  } view[] = {
    { 67.0, 133.0},
    { 90.0, 360.0},
    { 0.0, 0.0},
    {107.0, 277.0},
  };

  cout << "Ts =" << Ts << endl;
  cout << "dt =" << dt << endl;
  
  FILE *gp = popen("gnuplot","w");
  fprintf(gp, "set xrange[%lf:%lf]\n", -0.08, 0.10);
  fprintf(gp, "set yrange[%lf:%lf]\n", -0.02, 0.16);
  fprintf(gp, "set zrange[%lf:%lf]\n",  0.00, 0.20);
  //  fprintf(gp, "set view %lf, %lf\n"  , view[0].x, view[0].y);
  //  fprintf(gp, "set view %lf, %lf\n"  , view[1].x, view[1].y);
  fprintf(gp, "set view %lf, %lf\n"  , view[2].x, view[2].y);
  fprintf(gp, "set grid\n");

  //  fprintf(gp, "set terminal gif\n");

  for(int ii=10; ii < static_cast<int>(Ts/dt); ii+=10)
    { // Magick
      fprintf(gp, "splot ");
      fprintf(gp, "'dir-dat/%s_sinc_e1_%05d.dat' w l", prefix.c_str(), ii);
      fprintf(gp, ",");    
      fprintf(gp, "'dir-dat/%s_sinc_e2_%05d.dat' w l", prefix.c_str(), ii);
      fprintf(gp, "\n");    
      //    fprintf(gp, "pause 0.01\n");
      cout << __FUNCTION__<<"-"<< prefix <<"-"<< ii << endl;
    }

  for(int ii=10; ii < static_cast<int>(Ts/dt); ii+=10){
    //for(int ii=10; ii < static_cast<int>(3.0*ts[8]/dt); ii+=10){//Magick
    //fprintf(gp, "set output 'dir-gif/%s-%05d.gif\n", prefix.c_str(), ii);
    fprintf(gp, "splot ");
    fprintf(gp, "'dir-dat/%s_saiteki_e1_%05d.dat' w l", prefix.c_str(), ii);
    fprintf(gp, ",");
    fprintf(gp, "'dir-dat/%s_saiteki_e2_%05d.dat' w l ls 3", prefix.c_str(), ii);
    fprintf(gp, "\n");
    //    fprintf(gp, "pause 0.01\n");
    cout << __FUNCTION__<<"-"<< prefix <<"-"<< ii << endl;
  }
  pclose(gp);
  cout << __FUNCTION__<<"-"<< prefix << endl;
}

int
main(int argc, char *argv[])
{
  feenableexcept (FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  Ts = 0.64; t1 = Ts/8.0;  
  view_animation("curve");
  
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

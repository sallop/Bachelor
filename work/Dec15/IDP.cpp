#include "IDP.hpp"
#include "variables.hpp"
#include <cassert>
#include <cfloat>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Manipulator.hpp"
#include "IDP.hpp"

using namespace std;

void graph_volt_constituent(ofstream& ofs, double t, double ev,
			    double term1, double term2, double term3);
void graph_standard(ofstream& ofs, double t,
		    double z, double zv, double za,
		    double ev, double ia, double tau,
		    double energy, double ev_ia);

void graph_volt_constituent(FILE *fp, double t, double ev,
			    double term1, double term2, double term3);
void graph_standard(FILE *fp, double t,
		    double z, double zv, double za,
		    double ev, double ia, double tau,
		    double energy, double ev_ia);

void graph_endeffector(FILE *fp, VRManipulator &v1, Revolute &r3);



#define LEN(ary) (sizeof(ary) / sizeof(*ary) )
#define INITARY(ary) do{						\
    for(uint i=0; i < sizeof(ary)/sizeof(*ary); ++i){			\
      (ary)[i] = 0;							\
    }									\
  }while(0);

#undef DEBUG
#ifdef DEBUG
# define debug() cout << __FILE__<<"-"<<__FUNCTION__<<"-"<<__LINE__<< endl
#else
# define debug()
#endif
//static int IDP::arysz = 20;

IDP::IDP():v1(state_w.v1),
	   r1(state_w.r1), r2(state_w.r2), r3(state_w.r3),
	   e1(state_w.e1), e2(state_w.e2), e3(state_w.e3)// alias
{
  cout << __FUNCTION__ << endl;
}

IDP::IDP(string _name,
	 double _t_i, double _t_f,
	 double _dt, double _kiza,
	 int _n_repeat, int _n_divide, int _NN)
  : name(_name), t_i(_t_i), t_f(_t_f), dt(_dt), kiza(_kiza),
    n_repeat(_n_repeat), n_divide(_n_divide), NN(_NN),
    v1(state_w.v1),
    r1(state_w.r1), r2(state_w.r2), r3(state_w.r3),// alias
    e1(state_w.e1), e2(state_w.e2), e3(state_w.e3)// alias
{
  Ts = t_f - t_i;
  t1 = Ts/n_divide;
  for(uint i=0; i < LEN(ts); i++)
    {
      ts[i] = t1*i;
    }
  N3 = static_cast<int>(NN * 0.5 + 0.5);

  INITARY(Sth);
  INITARY(ntmax);
  INITARY(Qs);
  INITARY(Qv);
  INITARY(Qa);
  INITARY(at);
  //  INITARY(sat);
  for(int i = 0; i <= n_divide - 2; ++i)
    sat[i] = N3;

  for(int i = 0; i <= n_divide - 2; ++i)
    for(int j = 0; j <= NN; ++j)
      for(int k = 0; k <= NN; ++k)
	Mat[i][j][k] = 0;

  for(int i=0; i <= NN; ++i)
    for(int j=0; j <= NN; ++j)
      for(int k=0; k <= NN; ++k)
	Energ[i][j][k] = 0;

  
  //Enemin = new double(n_repeat); // hang in the air. offset
  for(int i=0; i < n_repeat; ++i)
    Enemin[i] = 0;
  
  // temp
  ntmax[0] = 0;
  ntmax[1] = static_cast<int>(ts[3]/dt);
  ntmax[2] = static_cast<int>(ts[4]/dt);
  ntmax[3] = static_cast<int>(ts[5]/dt);
  ntmax[4] = static_cast<int>(ts[8]/dt);
  // Qs_f, Qs_i => Init()

  sprintf(dname_dir_dat, "dir-dat");
  sprintf(dname_dir_etc, "dir-etc");
    
  sprintf(fname_sinc_st, "%s/%s_sinc_st.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_sinc_ev, "%s/%s_sinc_ev.dat"
	  , dname_dir_etc, name.c_str() );  
  sprintf(ftemp_saiteki_st, "%s/%s_saiteki_st_%%05d.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(ftemp_saiteki_ev, "%s/%s_saiteki_ev_%%05d.dat"
	  , dname_dir_etc, name.c_str() );  

  sprintf(fname_enemin2, "%s/%s_enemin2.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(ftemp_energy, "%s/%s_energy_%%05d.dat"
	  , dname_dir_etc, name.c_str() );

  sprintf(ftemp_sinc_e[1], "%s/%s_sinc_e1_%%05d.dat"
	  , dname_dir_dat, name.c_str() ); // 
  sprintf(ftemp_sinc_e[2], "%s/%s_sinc_e2_%%05d.dat"
	  , dname_dir_dat, name.c_str() ); // 
  sprintf(ftemp_sinc_e[3], "%s/%s_sinc_e3_%%05d.dat"
	  , dname_dir_dat, name.c_str() );
  // trajectory

  sprintf(ftemp_saiteki_e[1], "%s/%s_saiteki_e1_%%05d.dat"
	  , dname_dir_dat, name.c_str() );
  sprintf(ftemp_saiteki_e[2], "%s/%s_saiteki_e2_%%05d.dat"
	  , dname_dir_dat, name.c_str() );
  sprintf(ftemp_saiteki_e[3], "%s/%s_saiteki_e3_%%05d.dat"
	  , dname_dir_dat, name.c_str() );
  // trajectory

  for(uint i=0;
      i < sizeof(make_animation_flag)/sizeof(*make_animation_flag);
      ++i)
    {
      make_animation_flag[i] = 0;
    }
  for(uint i = 0; i < 10; i++)
    {
      make_animation_flag[i] = 1;
    }
  
  if (static_cast<int>(Ts/dt) < 100)
    {
      cout << "Sample points is less" << endl;
      assert(false);
    }

  printf("Ts/dt = %lf/%lf = %d\n", Ts, dt, static_cast<int>(Ts/dt) );
}

IDPUP::IDPUP(){}
IDPUP::~IDPUP(){}
IDPUP::IDPUP(string _name,
	     double _t_i, double _t_f,
	     double _dt, double _kiza,
	     int _n_repeat, int _n_divide, int _NN)
  :IDP( _name, _t_i,  _t_f, _dt,  _kiza, _n_repeat, _n_divide, _NN)
{
  debug();
}

IDPDOWN::IDPDOWN(){}
IDPDOWN::~IDPDOWN(){}
IDPDOWN::IDPDOWN(string _name,
		 double _t_i, double _t_f,
		 double _dt, double _kiza,
		 int _n_repeat, int _n_divide, int _NN)
  :IDP( _name, _t_i,  _t_f, _dt,  _kiza, _n_repeat, _n_divide, _NN)
{
}

IDPDIRECT::IDPDIRECT(){}
IDPDIRECT::~IDPDIRECT(){}
IDPDIRECT::IDPDIRECT(string _name,
		     double _t_i, double _t_f,
		     double _dt, double _kiza,
		     int _n_repeat, int _n_divide, int _NN)
  : IDP( _name,	 _t_i,  _t_f, _dt,  _kiza, _n_repeat, _n_divide, _NN)
{
  debug();
  sprintf(fname_r1_sinc_st, "%s/%s_r1_sinc_st.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r2_sinc_st, "%s/%s_r2_sinc_st.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r1_sinc_ev, "%s/%s_r1_sinc_ev.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r2_sinc_ev, "%s/%s_r2_sinc_ev.dat"
	  , dname_dir_etc, name.c_str() );

  sprintf(ftemp_r1_saiteki_st, "%s/%s_r1_saiteki_st_%%05d.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(ftemp_r2_saiteki_st, "%s/%s_r2_saiteki_st_%%05d.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(ftemp_r1_saiteki_ev, "%s/%s_r1_saiteki_ev_%%05d.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(ftemp_r2_saiteki_ev, "%s/%s_r2_saiteki_ev_%%05d.dat"
	  , dname_dir_etc, name.c_str() );

  sprintf(fname_r1_enemin2, "%s/%s_r1_enemin2.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r2_enemin2, "%s/%s_r2_enemin2.dat"
	  , dname_dir_etc, name.c_str() );
}

void IDPDIRECT::Init(//ScaraRobot _state_i,VRManipulator _v1_i,
		     Model _state_i,
		     SinCurv _sc,
		     Linear_curve_t _fk1,
		     Linear_curve_t _fk2,
		     double _Qs_i, double _Qs_f,
		     double _Qv_i, double _Qv_f,
		     double _Qa_i, double _Qa_f)  
{
  IDP::Init(_state_i,_sc,
	    _Qs_i, _Qs_f,
	    _Qv_i, _Qv_f,
	    _Qa_i, _Qa_f);
  fk1 = _fk1;
  fk2 = _fk2;
}

IDP::~IDP()
{
  //  delete Enemin;
}

void
IDP::Init(Model		_state_i,	//ScaraRobot _state_i,VRManipulator _v1_i,
	  SinCurv	_sc,
	  double	_Qs_i, double _Qs_f,
	  double	_Qv_i, double _Qv_f,
	  double	_Qa_i, double _Qa_f)
{
  //Qs_f   = ; => initialize at Inheritance
  // Qs_f ? Sth[n_repeat]
  state_i  = _state_i;
  state_w  = _state_i;		// now working states

  sc	   = _sc;
  Qs_i = _Qs_i;  Qs_f = _Qs_f;
  Qv_i = _Qv_i;  Qv_f = _Qv_f;
  Qa_i = _Qa_i;  Qa_f = _Qa_f;
  
  for(int i=0; i <= n_repeat - 2 ; ++i)
    {
      Sth[i] = sc.sin_th1(t1*i);
    }
}

void
IDP::start()
{
  debug();
  sinc();			// 初期軌道を与える
  ofstream ofs(fname_enemin2);
  for(l_nnn = 1; l_nnn <= n_repeat; ++l_nnn)
    {
      repeat();
      tansaku();
      //Enemin[l_nnn] = optimal_result(0);
      Enemin[l_nnn] = optimal_result( make_animation_flag[l_nnn] );

      ofs << l_nnn <<" "<< Enemin[l_nnn] <<" "<< kiza << endl;
      cout << name <<" "<<l_nnn <<" "<<Enemin[l_nnn]<<" "<< kiza << endl;
      if ((l_nnn > 15) && (Enemin[l_nnn] > Enemin[l_nnn - 1]))
	{
	  kiza *= 1.2;
	}
    }
  //state_f = state_w;
}

void
IDP::En_sim(int ntmax, int cur)  
{
  debug();
  double t = 0.0;
  double energy = 0.0;

  state_w = state_i;		// reinitialize
  for(int nt = 1; nt <= ntmax; nt++)
    {
      t += dt;
      sim(cur, t); //(cur, t)
      // inline masatu_jouken();
      // if (r3.zv < 0.0)
      // {
      // 	  e3.energy = DBL_MAX;
      // 	  break;
      // }
      energy = En();		// virtual function
    }
  Energ[at[cur] ][at[cur+1] ][at[cur+2] ] = energy;//=> sum_ev_ia_dt;
}

void IDP::zu(ostream& ofs)
{
  //cout << __FILE__ <<"-"<<__FUNCTION__<<"-"<< __LINE__ << endl;
  set_position();
  
  Circle3D real1_p0(r1.x0, r1.y0, r1.z0, 0.01);
  Circle3D real1_p1(r1.x , r1.y , r1.z0, 0.01);

  Circle3D real2_p0(r2.x0, r2.y0, r2.z0, 0.01);
  Circle3D real2_p1(r2.x , r2.y , r2.z0, 0.01);
  //  Circle3D real2_p2(r3.x0, r3.y0, r3.z0, 0.01);

  Circle3D real3_p0(r3.x0, r3.y0, r3.z0, 0.01);
  Circle3D real3_p1(r3.x0, r3.y0, r3.z , 0.01);
  
  Line3D real1_ln(r1.x0, r1.y0, r1.z0, r1.x,  r1.y , r1.z0 );
  //  Line3D real2_ln(r2.x0, r2.y0, r1.z0, r2.x,  r2.y , r1.z0 );
  Line3D real2_ln(r2.x0, r2.y0, r2.z0, r2.x,  r2.y , r2.z0 );
  Line3D real3_ln(r3.x0, r3.y0, r3.z0, r3.x0, r3.y0, r3.z  );
  
  Circle3D virt1_p0(v1.x0, v1.y0, v1.z0, 0.01);
  Circle3D virt1_p1(v1.x , v1.y , v1.z0, 0.01);

  Line3D virt1_ln(v1.x0, v1.y0, v1.z0,
		  v1.x , v1.y , v1.z0 );
  
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

void IDP::zu(ostream& ofs1, ostream& ofs2)
{
  //cout << __FILE__ <<"-"<<__FUNCTION__<<"-"<< __LINE__ << endl;
  set_position();
  
  Circle3D real1_p0(r1.x0, r1.y0, r1.z0, 0.01);
  Circle3D real1_p1(r1.x , r1.y , r1.z0, 0.01);

  Circle3D real2_p0(r2.x0, r2.y0, r2.z0, 0.01);
  Circle3D real2_p1(r2.x , r2.y , r2.z0, 0.01);
  //  Circle3D real2_p2(r3.x0, r3.y0, r3.z0, 0.01);

  Circle3D real3_p0(r3.x0, r3.y0, r3.z0, 0.01);
  Circle3D real3_p1(r3.x0, r3.y0, r3.z , 0.01);
  
  Line3D real1_ln(r1.x0, r1.y0, r1.z0, r1.x,  r1.y , r1.z0 );
  //  Line3D real2_ln(r2.x0, r2.y0, r1.z0, r2.x,  r2.y , r1.z0 );
  Line3D real2_ln(r2.x0, r2.y0, r2.z0, r2.x,  r2.y , r2.z0 );
  Line3D real3_ln(r3.x0, r3.y0, r3.z0, r3.x0, r3.y0, r3.z  );
  
  Circle3D virt1_p0(v1.x0, v1.y0, v1.z0, 0.01);
  Circle3D virt1_p1(v1.x , v1.y , v1.z0, 0.01);

  Line3D virt1_ln(v1.x0, v1.y0, v1.z0,
		  v1.x , v1.y , v1.z0 );
  
  // r1
 real1_p0.draw(ofs1);
 real1_p1.draw(ofs1);// real1_p2.draw(ofs);
 real1_ln.draw(ofs1);
  // r2
 real2_p0.draw(ofs1);
 real2_p1.draw(ofs1);// real2_p2.draw(ofs);
 real2_ln.draw(ofs1);
  // r3
 real3_p0.draw(ofs1);
 real3_p1.draw(ofs1);// real3_p2.draw(ofs);
 real3_ln.draw(ofs1);
  // v1
  virt1_p0.draw(ofs2);
  virt1_p1.draw(ofs2);// virt1_p2.draw(ofs);
  virt1_ln.draw(ofs2);  
}

void IDP::makeGifFile() const
{
  vector<string> fnames;
  for(uint i = 0; i < sizeof(make_animation_flag)/sizeof(*make_animation_flag); ++i)
    {
      if (make_animation_flag[i])
	{
	  char fname[32]; sprintf(fname, ftemp_energy, i);
	  fnames.push_back(fname);
	}
    }
  
  FILE* gp = popen("gnuplot", "w");
  fprintf(gp, "set terminal gif\n");
  fprintf(gp, "set grid\n");
  fprintf(gp, "set size 0.6, 0.6\n");
  fprintf(gp, "set output 'dir-gif/energy.gif'\n");

  fprintf(gp, "plot ");
  fprintf(gp, "'%s' u 1:2 w l ls 1 title 'sinc'", fnames[0].c_str() );
  for(uint i = 1; i < fnames.size() - 1; ++i)
    {
      fprintf(gp, ",'%s' u 1:2 w l ls 2 title '%d' ", fnames[i].c_str(), i);
    }
  fprintf(gp, ",'%s' u 1:2 w l ls 3 title 'result'\n", fnames.back().c_str() );  

  pclose(gp);
}

void IDPUP::sinc()
{
  char fname_energy[64]; sprintf(fname_energy, ftemp_energy, 0);
  ofstream ofs_st(fname_sinc_st), ofs_ev(fname_sinc_ev), ofs_energy(fname_energy);
  // FILE *fp_st = fopen(fname_sinc_st, "w");
  // FILE *fp_ev = fopen(fname_sinc_ev, "w");
  // FILE *fp_energy = fopen(fname_energy, "w");
  
  double t = 0;
  for(int i = 0; i <= static_cast<int>(Ts/dt); ++i)
    {
      graph_standard(ofs_st, t,
		     r3.z0, r3.zv, r3.za,
		     e3.ev, e3.ia, e3.tau,
		     e3.energy, e3.ev*e3.ia);

      graph_volt_constituent(ofs_ev, t, e3.ev,
			     b1*r3.Thv(), b2*r3.Tha(), b3*e3.tau);
      
      t = t + dt;
      r3.z0 = sc.sin_th1(t);
      r3.zv = sc.sin_th1v(t);
      r3.za = sc.sin_th1a(t);
      enemin2 = En();

      ofs_energy << t <<" "<< enemin2 << endl;
      //fprintf(ofs_energy, "%8.4lf %8.4lf\n", t, enemin2);
      
      if (i % 10 == 0)
	{
	  char fname[4][64];
	  sprintf(fname[1], ftemp_sinc_e[1], i);
	  sprintf(fname[2], ftemp_sinc_e[2], i);
	  // sprintf(fname[3], ftemp_sinc_e[3], i); !end-effector trajectory
	  ofstream ofs1(fname[1]);
	  ofstream ofs2(fname[2]);
	  // ofstream ofs3(fname[3]);
	  //zu( ofs1 );
	  zu(ofs1, ofs2);		// real, virtual
	  //zu(ofs1, ofs2, ofs3);
	}
      // output volt ampere field
    }
  enemin2 = 0.0;

  // ofs_st.close();
  //   ofs_ev.close();
  //ofs_energy.close();
  //  fclose(fp_st); fclose(fp_ev); fclose(fp_energy);
}

void IDPDOWN::sinc()
{
  char fname_energy[64]; sprintf(fname_energy, ftemp_energy, 0);
  ofstream ofs_st(fname_sinc_st), ofs_ev(fname_sinc_ev), ofs_energy(fname_energy);

  //  FILE	*fp2_st	   = fopen(fname_sinc_st,"w");
  //  FILE	*fp2_ev	   = fopen(fname_sinc_ev,"w");
  //FILE	*fp2_energy = fopen(fname_energy,"w");
  // FILE	*fp2_st	   =  fopen("dir-etc/fname_sinc_st","w");
  // FILE	*fp2_ev	   =  fopen("dir-etc/fname_sinc_ev","w");
  // FILE	*fp2_energy = fopen("dir-etc/fname_energy","w");
  
  double t = 0;  //  double s, v, a;
  for( int i = 0; i <= static_cast<int>(Ts/dt); ++i )
    {
      graph_standard(ofs_st, t,
		     r3.z0, r3.zv, r3.za,
		     e3.ev, e3.ia, e3.tau,
		     e3.energy, e3.ev*e3.ia);

      graph_volt_constituent(ofs_ev, t, e3.ev,
			     b1*r3.Thv(), b2*r3.Tha(), b3*e3.tau);
      
      t = t + dt;
      r3.z0 = sc.sin_th1(t);
      r3.zv = sc.sin_th1v(t);
      r3.za = sc.sin_th1a(t);
      enemin2 = En();

      ofs_energy << t <<" "<< enemin2 << endl;
      //      fprintf(fp2_energy, "%8.4lf %8.4lf\n", t, enemin2 );
      
      if (i % 10 == 0)
	{
	  char fname[4][64];
	  sprintf(fname[1], ftemp_sinc_e[1], i);
	  sprintf(fname[2], ftemp_sinc_e[2], i);
	  // sprintf(fname[3], ftemp_sinc_e[3], i);
	  // end-effector trajectory
	  ofstream ofs1(fname[1]);
	  ofstream ofs2(fname[2]);
	  // ofstream ofs3(fname[3]);
	  //zu( ofs1 ); all in
	  zu( ofs1, ofs2 );	// real, virt
	  //zu( ofs1, ofs2, ofs3 );// real, virt, trjc
	}
      // output volt ampere field
      // plot("volt ampere field");
    }
  //  fclose(fp2_st);
  //   fclose(fp2_ev);
  // fclose(fp2_energy);  
}

void IDPDIRECT::sinc()
{
  char fname_energy[64]; sprintf(fname_energy, ftemp_energy, 0);
  ofstream ofs_energy(fname_energy);
  ofstream ofs_r1_st(fname_r1_sinc_st), ofs_r1_ev(fname_r1_sinc_ev);
  ofstream ofs_r2_st(fname_r2_sinc_st), ofs_r2_ev(fname_r2_sinc_ev);
  
  //  FILE* fp_energy = fopen(fname_energy,"w");
  printf("line=%d %s\n", __LINE__, fname_r1_sinc_st);
  printf("line=%d %s\n", __LINE__, fname_r2_sinc_st);
  printf("line=%d %s\n", __LINE__, fname_r1_sinc_ev);
  printf("line=%d %s\n", __LINE__, fname_r2_sinc_ev);

  //   FILE* fp_r1_st = fopen(fname_r1_sinc_st,"w");
  //   FILE* fp_r1_ev = fopen(fname_r1_sinc_ev,"w");
  //   FILE* fp_r2_st = fopen(fname_r2_sinc_st,"w");
  //   FILE* fp_r2_ev = fopen(fname_r2_sinc_ev,"w");
  
  double t = 0;
  for( int i = 0; i <= static_cast<int>(Ts/dt); ++i )
    {
      graph_standard(ofs_r1_st, t,
		     r1.th, r1.thv, r1.tha,
		     e1.ev, e1.ia, e1.tau,
		     e1.energy, e1.ev*e1.ia);
      graph_volt_constituent(ofs_r1_ev, t, e1.ev,
			     b1*r1.thv, b2*r1.tha, b3*e1.tau);
      
      graph_standard(ofs_r2_st, t,
		     r2.th, r2.thv, r2.tha,
		     e2.ev, e2.ia, e2.tau,
		     e2.energy, e2.ev*e2.ia);
      
      graph_volt_constituent(ofs_r2_ev, t, e2.ev,
			     b1*r2.thv, b2*r2.tha, b3*e2.tau);

      t = t + dt;
      v1.th  = sc.sin_th1(t);
      v1.thv = sc.sin_th1v(t);
      v1.tha = sc.sin_th1a(t);
      enemin2 = En();

      ofs_energy << t <<" "<< enemin2 << endl;
      //      fprintf(fp_energy, "%8.4lf %8.4lf\n", t, enemin2 );

      if (i % 10 == 0)
	{
	  char fname[4][64];
	  sprintf(fname[1], ftemp_sinc_e[1], i);
	  sprintf(fname[2], ftemp_sinc_e[2], i);	
	  // sprintf(fname[3], fname_sinc_e[3],  name.c_str(), i);
	  ofstream ofs1(fname[1]);
	  ofstream ofs2(fname[2]);
	  // ofstream ofs3(fname[3]);
	  //zu( ofs1 ); all in
	  zu( ofs1, ofs2 );
	  //	zu( ofs1, ofs2, ofs3 ); !not implement
	}
      // output volt ampere field
      // plot("");
    }
  
  //   fclose(fp_energy);
  //   fclose(fp_r1_st);  fclose(fp_r1_ev);
  //   fclose(fp_r2_st);  fclose(fp_r2_ev);
}

void IDP::repeat()		// ofstream
{
  int ff1;
  //Magick
  for(int ii=6; ii > 0; --ii)
    {  /* 中心位置に探索結果の最適経路を選ぶ */
      Sth[ii] = Sth[ii] + kiza*(sat[ii] - N3);
    }

  {// get summation
    double sum = 0.0;
    for(int ii=6; ii > 0; --ii){//Magick
      sum = sum + fabs(sat[ii] - N3);
    }
    ff1 = static_cast<int>(sum);
  }
  
  if (0 == ff1){
    kiza = kiza*0.8;
    printf("%s-%s, kiza1=%lf, nnn=%d\n", __FILE__, __FUNCTION__, kiza, l_nnn);
  }
}

void IDP::tansaku()
{
  // at[cur] = 1 : NN;
  // cur = 1 : (n_divide - 2) - 2;
  for(int cur = 1; cur <= 4; cur++){//Magick latter, replace the divide number
    for(at[cur+2] = 1; at[cur+2] <= NN; at[cur+2]++){
      for(at[cur+1]=1; at[cur+1] <= NN; at[cur+1]++){
	for(at[cur]=1; at[cur  ] <= NN; at[cur  ]++){

	  Qs[cur+2] = Sth[cur+2] + kiza*(at[cur+2] - N3);	  
	  Qs[cur+1] = Sth[cur+1] + kiza*(at[cur+1] - N3);      	  
	  Qs[cur  ] = Sth[cur  ] + kiza*(at[cur  ] - N3);

	  for(int i = cur; --i > 0;){
	    at[i]  = Mat[i][at[i+1] ][at[i+2] ];
	    Qs[i] = Sth[i] + kiza*(at[i] - N3);
	  }
	  settei(cur);// set  Qv[], Qa[];
	  En_sim( ntmax[cur], cur );
	}
      }
    }
    Hikaku(cur);
  }
  return;
}

double
IDP::optimal_result(int make_animation_flag)
{
  int cur;
  double t = 0;

  // last 3 step(6,5,4) is conserved at Hikaku1.
  //  for(int ii = 6; ii > 3; ii--){ at[ii] = sat[ii]; }//Magick  
  
  //for(int ii = (n_divide -2); ii > ((n_divide-2) - 3); --ii)
  //  {
  //    at[ii] = sat[ii];
  //  }
  at[6] = sat[6];
  at[5] = sat[5];
  at[4] = sat[4];
  
  //for(int ii = 3; ii > 0; ii--)
  for(int ii = ((n_divide-2)-3); ii > 0; --ii){//Magick
    at[ii] = Mat[ii][at[ii+1] ][at[ii+2] ];
    sat[ii] =  at[ii];
  }  
  //  for(int ii = 6; ii > 0; ii--){//Magick
  for(int ii = n_divide-2; ii > 0; --ii){//Magick
      Qs[ii] = Sth[ii] + kiza*(at[ii] - N3);
  }

  //cur = 4; ((8 - 2) - 3) + 1//n_divide=8
  cur = ((n_divide - 2) - 3) + 1;
  settei(cur);
  //    reinitialize();
  ofstream ofs_st, ofs_ev, ofs_energy;
  if (make_animation_flag)
    {
      char fname[64];
      sprintf(fname, ftemp_saiteki_st, l_nnn); ofs_st.open(fname);
      sprintf(fname, ftemp_saiteki_ev, l_nnn); ofs_ev.open(fname);
      sprintf(fname, ftemp_energy, l_nnn); ofs_energy.open(fname);
    }

  //  for(int nt = 1; nt <= ntmax[ cur ]; nt++)
  for(int nt = 1; nt <= ntmax[ cur ]; nt++)
    {  
      t += dt;
      sim(cur, t);
      enemin2 = En();
      
      ofs_energy << t <<" "<<enemin2<< endl;
      
      if (make_animation_flag)
	{
	  graph_standard(ofs_st, t,
			 r3.z, r3.zv, r3.za,
			 e3.ev, e3.ia, e3.tau,
			 e3.energy, e3.ev*e3.ia);
      
	  graph_volt_constituent(ofs_ev, t, e3.ev,
				 b1*r3.Thv(),
				 b2*r3.Tha(),
				 b3*e3.tau);
	}// fi make_animation_flag

      if (nt % 10 == 0)
	{
	  char fname[4][64];
	  sprintf(fname[1], ftemp_saiteki_e[1], nt);
	  sprintf(fname[2], ftemp_saiteki_e[2], nt);
	  sprintf(fname[3], ftemp_saiteki_e[3], nt);
	  ofstream ofs1(fname[1]);
	  ofstream ofs2(fname[2]);
	  //	  ofstream ofs3(fname[3]);
	  zu(ofs1, ofs2);
	}
    } //endfor

  if (make_animation_flag)
    {
      ofs_st.close();
      ofs_ev.close();
      ofs_energy.close();
    }
  
  return enemin2;
}

void IDP::settei(int cur)
{
  //  assert(cur <= (n_divide-2)-2);
  Qv[1] = 2.0*(Qs[1] - Qs_i)/t1;
  Qa[1] = Qv[1]/t1;
  Qv[2] = 2.0*(Qs[2] - Qs[1])/t1 - Qv[1];
  Qa[2] = (Qv[2] - Qv[1])/t1;
  Qv[3] = 2.0*(Qs[3] - Qs[2])/t1 - Qv[2];
  Qa[3] = (Qv[3] - Qv[2])/t1;
  if (cur == 1) return;
    
  Qv[4] = 2.0*(Qs[4] - Qs[3])/t1 - Qv[3];
  Qa[4] = (Qv[4] - Qv[3])/t1;
  if (cur == 2) return;
    
  Qv[5] = 2.0*(Qs[5] - Qs[4])/t1 - Qv[4];
  Qa[5] = (Qv[5] - Qv[4])/t1;
  if (cur == 3) return;
    
  Qv[6] = 2.0*(Qs[6] - Qs[5])/t1 - Qv[5];
  Qa[6] = (Qv[6] - Qv[5])/t1;
  // will become if (step == 4) return;
  
  Qs[8] = Qs_f;
  Qv[8] = 0.0; // vary with condition
  Qv[7] = (Qs[8] - Qs[6])/t1 - (Qv[6]+Qv[8])/2.0;
  Qa[7] = (Qv[7] - Qv[6])/t1;
  Qs[7] = Qs[6] + Qv[6] * t1 + Qa[7]*t1*t1/2.0;
  Qa[8] = (Qv[8] - Qv[7])/t1;
  if (cur == 4) return;

  // if (cur == 5?)
  //  return; ((n_divide-2)-3) + 1 < mat_max
  // if (cur <= mat_max)
  //     settei1()
  // else
  //     settei2(get the last two step acceration and velocity)
  // 
  // mat_max = ((n_divide - 2) - 3) + 1 (n_divide=8; 4);
}

// inline
//void IDP::sim(int cur_step, double& refs, double& refv, double& refa)
//{
//  refa = Qa[cur_step];
//  refv = Qv[cur_step - 1]
//    + refa*(t - ts[cur_step - 1]);
//  refs = Qs[cur_step - 1]
//    + Qv[cur_step - 1]*(t - ts[cur_step - 1])
//    + 0.5*refa*SQ(t - ts[cur_step - 1]);
//}

void IDP::sim(int cur, double t)//Magick_numbers
{
  //  assert(cur <= (n_divide-2)-2);
#define SQ(x) ((x)*(x))
  if (t < ts[1])
    {
      r3.za = Qa[1];
      r3.zv = r3.za*t;
      r3.z0 = 0.5*r3.za*t*t + Qs_i;
    }
  if (ts[1] < t && t < ts[2])
    {
      r3.za = Qa[2];
      r3.zv = Qv[1] + r3.za*(t - ts[1]);
      r3.z0 = Qs[1] + Qv[1]*(t - ts[1]) + 0.5*r3.za*SQ(t - ts[1]);
      //r3.z0 = th1[1] + th1v[1]*(t - ts[1]) + 0.5*r3.za*pow((t - ts[1]), 2);	
    }
  if (ts[2] < t && t < ts[3])
    {
      r3.za = Qa[3];
      r3.zv = Qv[2] + r3.za*(t - ts[2]);
      //	r3.z0 = th1[2] + th1v[2]*(t - ts[2]) + 0.5*r3.za*pow((t - ts[2]), 2);
      r3.z0 = Qs[2] + Qv[2]*(t - ts[2]) + 0.5*r3.za*SQ(t - ts[2]);	
    }
  if (1 == cur)
    return;

  if (ts[3] < t && t < ts[4])
    {
      r3.za = Qa[4];
      r3.zv = Qv[3] + r3.za*(t - ts[3]);
      r3.z0 = Qs[3] + Qv[3]*(t - ts[3]) + 0.5*r3.za*SQ(t - ts[3]);
      //r3.z0 = th1[3] + th1v[3]*(t - ts[3]) + 0.5*r3.za*pow((t - ts[3]), 2);
    }
  if (2 == cur)
    return;

  if (ts[4] < t && t < ts[5])
    {
      r3.za = Qa[5];
      r3.zv = Qv[4] + r3.za*(t - ts[4]);
      r3.z0  = Qs[4] + Qv[4]*(t - ts[4]) + 0.5*r3.za*SQ(t - ts[4]);
      //	r3.z0  = th1[4] + th1v[4]*(t - ts[4]) + 0.5*r3.za*pow((t - ts[4]), 2);
    }
  if (3 == cur)
    return;

  if (ts[5] < t && t < ts[6])
    {
      r3.za = Qa[6];
      r3.zv = Qv[5] + r3.za*(t - ts[5]);
      r3.z0 = Qs[5] + Qv[5]*(t - ts[5]) + 0.5*r3.za*SQ(t - ts[5]);
      //r3.z0 = th1[5] + th1v[5]*(t - ts[5]) + 0.5*r3.za*pow((t - ts[5]), 2);
    }
  if (ts[6] < t && t < ts[7])
    {
      r3.za = Qa[7];
      r3.zv = Qv[6] + r3.za*(t - ts[6]);
      r3.z0 = Qs[6] + Qv[6]*(t - ts[6]) + 0.5*r3.za*SQ(t - ts[6]);
      //	r3.z0 = th1[6] + th1v[6]*(t - ts[6]) + 0.5*r3.za*pow((t - ts[6]), 2);
    }
  if (ts[7] < t && t < ts[8])
    {
      r3.za = Qa[8];
      r3.zv = Qv[7] + r3.za*(t - ts[7]);
      r3.z0 = Qs[7] + Qv[7]*(t - ts[7]) + 0.5*r3.za*SQ(t - ts[7]);
      //	r3.z0  = th1[7] + th1v[7]*(t - ts[7]) + 0.5*r3.za*pow((t - ts[7]), 2);
    }
  if (ts[8] <= t)
    {// Through path
    }
  if (4 == cur)
    return;
#undef SQ
}

void IDPUP::sim(int cur, double t)//Magick_numbers
{
  //  assert(cur <= (n_divide-2)-2);
#define SQ(x) ((x)*(x))
  if (t < ts[1])
    {
      r3.za = Qa[1];
      r3.zv = r3.za*t;
      r3.z0 = 0.5*r3.za*t*t + Qs_i;
    }
  if (ts[1] < t && t < ts[2])
    {
      r3.za = Qa[2];
      r3.zv = Qv[1] + r3.za*(t - ts[1]);
      r3.z0 = Qs[1] + Qv[1]*(t - ts[1]) + 0.5*r3.za*SQ(t - ts[1]);
    }
  if (ts[2] < t && t < ts[3])
    {
      r3.za = Qa[3];
      r3.zv = Qv[2] + r3.za*(t - ts[2]);
      r3.z0 = Qs[2] + Qv[2]*(t - ts[2]) + 0.5*r3.za*SQ(t - ts[2]);	
    }
  if (1 == cur)
    return;

  if (ts[3] < t && t < ts[4])
    {
      r3.za = Qa[4];
      r3.zv = Qv[3] + r3.za*(t - ts[3]);
      r3.z0 = Qs[3] + Qv[3]*(t - ts[3]) + 0.5*r3.za*SQ(t - ts[3]);
    }
  if (2 == cur)
    return;

  if (ts[4] < t && t < ts[5])
    {
      r3.za = Qa[5];
      r3.zv = Qv[4] + r3.za*(t - ts[4]);
      r3.z0  = Qs[4] + Qv[4]*(t - ts[4]) + 0.5*r3.za*SQ(t - ts[4]);
    }
  if (3 == cur)
    return;

  if (ts[5] < t && t < ts[6])
    {
      r3.za = Qa[6];
      r3.zv = Qv[5] + r3.za*(t - ts[5]);
      r3.z0 = Qs[5] + Qv[5]*(t - ts[5]) + 0.5*r3.za*SQ(t - ts[5]);
    }
  if (ts[6] < t && t < ts[7])
    {
      r3.za = Qa[7];
      r3.zv = Qv[6] + r3.za*(t - ts[6]);
      r3.z0 = Qs[6] + Qv[6]*(t - ts[6]) + 0.5*r3.za*SQ(t - ts[6]);
    }
  if (ts[7] < t && t < ts[8])
    {
      r3.za = Qa[8];
      r3.zv = Qv[7] + r3.za*(t - ts[7]);
      r3.z0 = Qs[7] + Qv[7]*(t - ts[7]) + 0.5*r3.za*SQ(t - ts[7]);
    }
  if (ts[8] <= t)
    {// Through path
    }
  if (4 == cur)
    return;
#undef SQ
}

void IDPDOWN::sim(int cur, double t)//Magick_numbers
{
  //  assert(cur <= (n_divide-2)-2);
#define SQ(x) ((x)*(x))
  if (t < ts[1])
    {
      r3.za = Qa[1];
      r3.zv = r3.za*t;
      r3.z0 = 0.5*r3.za*t*t + Qs_i;
    }
  if (ts[1] < t && t < ts[2])
    {
      r3.za = Qa[2];
      r3.zv = Qv[1] + r3.za*(t - ts[1]);
      r3.z0 = Qs[1] + Qv[1]*(t - ts[1]) + 0.5*r3.za*SQ(t - ts[1]);
    }
  if (ts[2] < t && t < ts[3])
    {
      r3.za = Qa[3];
      r3.zv = Qv[2] + r3.za*(t - ts[2]);
      r3.z0 = Qs[2] + Qv[2]*(t - ts[2]) + 0.5*r3.za*SQ(t - ts[2]);	
    }
  if (1 == cur)
    return;

  if (ts[3] < t && t < ts[4])
    {
      r3.za = Qa[4];
      r3.zv = Qv[3] + r3.za*(t - ts[3]);
      r3.z0 = Qs[3] + Qv[3]*(t - ts[3]) + 0.5*r3.za*SQ(t - ts[3]);
    }
  if (2 == cur)
    return;

  if (ts[4] < t && t < ts[5])
    {
      r3.za = Qa[5];
      r3.zv = Qv[4] + r3.za*(t - ts[4]);
      r3.z0  = Qs[4] + Qv[4]*(t - ts[4]) + 0.5*r3.za*SQ(t - ts[4]);
    }
  if (3 == cur)
    return;

  if (ts[5] < t && t < ts[6])
    {
      r3.za = Qa[6];
      r3.zv = Qv[5] + r3.za*(t - ts[5]);
      r3.z0 = Qs[5] + Qv[5]*(t - ts[5]) + 0.5*r3.za*SQ(t - ts[5]);
    }
  if (ts[6] < t && t < ts[7])
    {
      r3.za = Qa[7];
      r3.zv = Qv[6] + r3.za*(t - ts[6]);
      r3.z0 = Qs[6] + Qv[6]*(t - ts[6]) + 0.5*r3.za*SQ(t - ts[6]);
    }
  if (ts[7] < t && t < ts[8])
    {
      r3.za = Qa[8];
      r3.zv = Qv[7] + r3.za*(t - ts[7]);
      r3.z0 = Qs[7] + Qv[7]*(t - ts[7]) + 0.5*r3.za*SQ(t - ts[7]);
    }
  if (ts[8] <= t)
    {// Through path
    }
  if (4 == cur)
    return;
#undef SQ
}

// 上のモジュールでtの値を判定して、添字のインデックスを渡す
// ここのモジュールでは、
// Qa[step];
// Qv[step]+v1.tha*(t - ts[step-1]);
// Qs[step-1]+Qv[step-1]*(t - ts[step-1]) + 0.5*v1.tha*SQ(t - ts[step - 1];
// だけをもとめ、呼び出し側で引数にstepを指定するように変更する。

void IDPDIRECT::sim(int cur, double t)//Magick_numbers
{
  //  assert(cur <= (n_divide-2)-2);
#define SQ(x) ((x)*(x))
  if (t < ts[1])
    {
      v1.tha = Qa[1];
      v1.thv = v1.tha*t;
      v1.th = 0.5*v1.tha*t*t + Qs_i;
    }
  if (ts[1] < t && t < ts[2])
    {
      v1.tha = Qa[2];
      v1.thv = Qv[1] + v1.tha*(t - ts[1]);
      v1.th = Qs[1] + Qv[1]*(t - ts[1]) + 0.5*v1.tha*SQ(t - ts[1]);
    }
  if (ts[2] < t && t < ts[3])
    {
      v1.tha = Qa[3];
      v1.thv = Qv[2] + v1.tha*(t - ts[2]);
      v1.th = Qs[2] + Qv[2]*(t - ts[2]) + 0.5*v1.tha*SQ(t - ts[2]);	
    }
  if (1 == cur)
    return;

  if (ts[3] < t && t < ts[4])
    {
      v1.tha = Qa[4];
      v1.thv = Qv[3] + v1.tha*(t - ts[3]);
      v1.th = Qs[3] + Qv[3]*(t - ts[3]) + 0.5*v1.tha*SQ(t - ts[3]);
    }
  if (2 == cur)
    return;

  if (ts[4] < t && t < ts[5])
    {
      v1.tha = Qa[5];
      v1.thv = Qv[4] + v1.tha*(t - ts[4]);
      v1.th  = Qs[4] + Qv[4]*(t - ts[4]) + 0.5*v1.tha*SQ(t - ts[4]);
    }
  if (3 == cur)
    return;

  if (ts[5] < t && t < ts[6])
    {
      v1.tha = Qa[6];
      v1.thv = Qv[5] + v1.tha*(t - ts[5]);
      v1.th = Qs[5] + Qv[5]*(t - ts[5]) + 0.5*v1.tha*SQ(t - ts[5]);
    }
  if (ts[6] < t && t < ts[7])
    {
      v1.tha = Qa[7];
      v1.thv = Qv[6] + v1.tha*(t - ts[6]);
      v1.th = Qs[6] + Qv[6]*(t - ts[6]) + 0.5*v1.tha*SQ(t - ts[6]);
    }
  if (ts[7] < t && t < ts[8])
    {
      v1.tha = Qa[8];
      v1.thv = Qv[7] + v1.tha*(t - ts[7]);
      v1.th = Qs[7] + Qv[7]*(t - ts[7]) + 0.5*v1.tha*SQ(t - ts[7]);
    }
  if (ts[8] <= t)
    {// Through path
    }
  if (4 == cur)
    return;
#undef SQ
}

double IDP::Hikaku(int current_step)
{
  //assert(current_step <= (n_divide-2)-2);
  int na1, na2, na3;//, na4;
  enemin2 = DBL_MAX;
  for(na3=1; na3<=NN; na3++)
    {
      for(na2=1; na2<=NN; na2++)
	{
	  enemin1 = DBL_MAX;
	  /* -------------------- */
	  for(na1=1; na1<=NN; na1++)
	    {
	      if (enemin1 > Energ[na1][na2][na3])
		{
		  enemin1 = Energ[na1][na2][na3];
		  Mat[current_step][na2][na3] = na1;
		}
	
	      if (enemin2 > Energ[na1][na2][na3])
		{
		  enemin2 = Energ[na1][na2][na3];
		  sat[6] = na3; //Magick max_cur(= n_divide-2-2)
		  sat[5] = na2; //Magick max_cur-1
		  sat[4] = na1; //Magick max_cur
		  //sat[n_divide - 2] = na3; //Magick max_cur(= n_divide-2-2)
		  //sat[n_divide - 1] = na2; //Magick max_cur-1
		  //sat[n_divide - 0] = na1; //Magick max_cur
		}
	    }
	}
    }
  //assert(enemin1 >= 0 && enemin2 >= 0);
  return enemin2;
}

#define DISPLAY(out, num) (out) << #num "=" << (num) << endl
#define ARYDISPLAY(out, ary) do{			\
    for(int ii = 0; ii < arysz; ++ii){			\
      (out) << #ary "["<<ii<<"]=" << (ary)[ii] << endl;	\
    }							\
  } while(0)

void IDP::dump(ostream& os)
{
  DISPLAY(os, n_divide);
  DISPLAY(os, dt);
  DISPLAY(os, Ts);
  DISPLAY(os, t1);
  ARYDISPLAY(os, ts);// SIGFPE
  DISPLAY(os, Qs_i);	// not initialize
  DISPLAY(os, Qs_f);
  DISPLAY(os, Qv_i);	// not initialize
  DISPLAY(os, Qv_f);
  DISPLAY(os, Qa_i);	// not initialize
  DISPLAY(os, Qa_f);    
  // Program terminated with SIGFPE,
  // Arithmetic exception.
  ARYDISPLAY(os, Sth );
  //ARYDISPLAY(os, th );
  ARYDISPLAY(os, ntmax );
  DISPLAY(os, NN);
  DISPLAY(os, N3);
  ARYDISPLAY(os, Qs);
  ARYDISPLAY(os, Qv);
  ARYDISPLAY(os, Qa);    
  DISPLAY(os, kiza);
  //DISPLAY(os, t);
}
#undef DISPLAY
#undef ARYDISPLAY

const Model& IDP::getState() const
{
  return state_w;
}


double IDPUP::En()
{
  e3.tau = r3.m*(r3.za + g)*r3.radius;// pulley's radius and torque.
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
  
  if (e3.ev*e3.ia > 0.0){
    e3.energy = e3.energy + e3.ev*e3.ia*dt;
  }
  return e3.energy;
}

double IDPDOWN::En()
{
  e3.tau = r3.m*(r3.za + g)*r3.radius;// pulley's radius and torque.
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
  return e3.energy;
}

double IDPDIRECT::det(double a11, double a12, double a21, double a22)
{
  return a11*a22 - a12*a21;
}

double IDPDIRECT::En()
{
  double detJ, J11, J12, J21, J22, coef;
  double cth1, cth12, cth2;
  double sth1, sth12, sth2;  
  double keirosokudo, PL3, XXa, YYa;
  //  double fxx, fyy;

  cth1 = cos(r1.th); cth12 = cos(r1.th + r2.th); cth2 = cos(r2.th);
  sth1 = sin(r1.th); sth12 = sin(r1.th + r2.th); sth2 = sin(r2.th);
  
  fk1.a = (pnt_f.y - pnt_i.y)/(pnt_f.x - pnt_i.x);// constant
  fk1.b = pnt_i.y;				  // constant
  fk2.a = tan(v1.th);		// fk2.a = tan(th_.s); variable
  fk2.b = v1.y0 - fk2.a*v1.x0 ; // => check x0         variable
  v1.x  = (fk2.b - fk1.b)/(fk1.a - fk2.a);
  v1.y  = fk2.a*v1.x + fk2.b;

  v1.l = hypot(v1.x - v1.x0 , v1.y - v1.y0);
  coef = ( fk1.a*tan(v1.th) + 1.0 )/( fk1.a - tan(v1.th) );
  v1.lv = coef*v1.l*v1.thv;  // v1.lv = coef*v1.l*th_.v;
  v1.xv = v1.lv*cos(v1.th) - v1.l*v1.thv*sin(v1.th);
  v1.yv = v1.lv*sin(v1.th) + v1.l*v1.thv*cos(v1.th);

  keirosokudo = sqrt(v1.xv*v1.xv + v1.yv*v1.yv);
  J11 = -r1.l*sth1 - r2.l*sth12;
  J12 = -r2.l*sth12;
  J21 =  r1.l*cth1 + r2.l*cth12;
  J22 =  r2.l*cth12;

  detJ = det(J11, J12, J21, J22);
  r1.thv = det(v1.xv, J12  , v1.yv, J22  )/detJ;
  // => arithmetic exception
  r2.thv = det(J11  , v1.xv, J21  , v1.yv)/detJ;

  PL3 = sqrt(v1.x*v1.x + v1.y*v1.y); // v1.x ? (v1.x - v1.x0)
  r1.th = atan2(v1.y, v1.x) - acos(0.5*PL3/r1.l);
  r2.th = 2.0*acos(0.5*PL3/r1.l);


  v1.la = coef*(v1.l*v1.tha + 2.0*v1.lv*v1.thv) + v1.l*v1.thv*v1.thv;
  v1.xa  = v1.la*cos(v1.th) - 2.0*v1.lv*v1.thv*sin(v1.th)
    - v1.l*v1.tha*sin(v1.th) - v1.l*v1.thv*v1.thv*cos(v1.th);


  v1.ya  = v1.la*sin(v1.th) + 2.0*v1.lv*v1.thv*cos(v1.th)
    + v1.l*v1.tha*cos(v1.th) - v1.l*v1.thv*v1.thv*sin(v1.th);
  
  XXa = v1.xa + r1.thv*r1.thv*r1.l*cth1
    + (r1.thv + r2.thv)*(r1.thv + r2.thv)*r2.l*cth12;
  YYa = v1.ya + r1.thv*r1.thv*r1.l*sth1
    + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*sth12;

  r1.tha = det(XXa, J12, YYa, J22)/detJ;
  r2.tha = det(J11, XXa, J21, YYa)/detJ;

  cth1 = cos(r1.th); cth12 = cos(r1.th + r2.th); cth2 = cos(r2.th);
  sth1 = sin(r1.th); sth12 = sin(r1.th + r2.th); sth2 = sin(r2.th);
  
  e1.tau = (a1+a2+a3+2.0*a4*cth2)*r1.tha + (a2+a3+a4*cth2)*r2.tha 
    - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*sth2;
  e2.tau = (a2 + a3 + a4*cth2)*r1.tha
    + (a2 + a3)*r2.tha
    + a4*r1.thv*r1.thv*sth2;

  e1.ev = b1*r1.thv + b2*r1.tha + b3*e1.tau;
  e1.ia = (e1.ev - Kv*r1.thv)/Ra;
  e2.ev = b1*r2.thv + b2*r2.tha + b3*e2.tau;
  e2.ia = (e2.ev - Kv*r2.thv)/Ra;
    
  //if (ev1*ia1 > 0.0){ Ene_1 = Ene_1 + ev1*ia1*dt; }
  //if (ev2*ia2 > 0.0){ Ene_2 = Ene_2 + ev2*ia2*dt; }
#define WATT( e )  ( (e).ev*(e).ia )
  if ( WATT(e1) > 0.0){
    e1.energy = e1.energy + WATT(e1)*dt;
  }
  if ( WATT(e2) > 0.0){ //      Ene_2 = Ene_2 + ev2*ia2*dt;
    e2.energy = e2.energy + WATT(e2)*dt;
  }
#undef WATT
  //Ene123 = e1.energy + e2.energy + e3.energy;
  //  return Ene123; // => Nov 23
  //Ene123 = e1.energy + e2.energy;
  return e1.energy + e2.energy;
}

void
IDP::set_position(){}

void
IDPUP::set_position()
{
  //  v1.x0 = pnt_v0.x; v1.y0 = pnt_v0.y;
  v1.z0 = r3.z0  ;
  //  v1.x  = pnt_i.x; v1.y  = pnt_i.y;
  v1.x = v1.x0 + v1.l*cos(v1.th);
  v1.y = v1.y0 + v1.l*sin(v1.th); 
  r3.z = r3.z0 + r3.l;  
}

void
IDPDOWN::set_position()
{
  //  v1.x0 = pnt_v0.x; v1.y0 = pnt_v0.y;
  v1.z0 = r3.z0  ;
  //  v1.x  = pnt_i.x; v1.y  = pnt_i.y;
  v1.x = v1.x0 + v1.l*cos(v1.th);
  v1.y = v1.y0 + v1.l*sin(v1.th); 
  r3.z = r3.z0 + r3.l;  
}

void
IDPDIRECT::set_position()
{
  r1.x = r1.x0 + r1.l*cos(r1.th);
  r1.y = r1.y0 + r1.l*sin(r1.th);
  //  r1.z = r1.z0 ;
  r2.x0 = r1.x ; r2.x = r2.x0 + r2.l*cos( r1.th + r2.th );
  r2.y0 = r1.y ; r2.y = r2.y0 + r2.l*sin( r1.th + r2.th );
  r2.z0 = r1.z0;

  r3.x0 = r2.x ; r3.y0 = r2.y;

  v1.x = v1.x0 + v1.l*cos(v1.th);
  v1.y = v1.y0 + v1.l*sin(v1.th);
  //r3.z0 = r2.z0; r3.z  = r3.z0 + r3.l;
}

Model::Model(){}
Model::Model(VRManipulator _v1,
	     Revolute _r1, Revolute _r2, Prismatic _r3,
	     Energy _e1, Energy _e2, Energy _e3)
  : v1(_v1),
    r1(_r1), r2(_r2), r3(_r3),
    e1(_e1), e2(_e2), e3(_e3)
{
}
Model::~Model(){}

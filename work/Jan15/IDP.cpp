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
			    double term1 , double term2, double term3);

void graph_standard(ofstream& ofs, double t    ,
		    double z     , double zv   , double za ,
		    double ev    , double ia   , double tau,
		    double energy, double ev_ia);

void graph_endeffector(ofstream& ofs,
		       double t, VRManipulator &v1, Prismatic &r3);
void graph_endeffector(FILE *fp, double t, VRManipulator &v1, Prismatic &r3);

void graph_volt_constituent(FILE *fp, double t, double ev,
			    double term1, double term2, double term3);

void graph_standard(FILE *fp, double t ,
		    double z     , double zv, double za ,
		    double ev    , double ia, double tau,
		    double energy, double ev_ia);

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
    e1(state_w.e1), e2(state_w.e2), e3(state_w.e3) // alias
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
    {
      sat[i] = N3;
    }
  
  for(int i = 0; i <= n_divide - 2; ++i)
    for(int j = 0; j <= NN; ++j)
      for(int k = 0; k <= NN; ++k)
	{
	  Mat[i][j][k] = 0;
	}
  for(int i=0; i <= NN; ++i)
    for(int j=0; j <= NN; ++j)
      for(int k=0; k <= NN; ++k)
	{
	  Energ[i][j][k] = 0;
	}
  //Enemin = new double(n_repeat); // hang in the air. offset
  //for(int i=0; i < n_repeat; ++i)
  //{
  //Enemin[i] = 0;
  //  }
  
  for(int j = 0; j < 4; ++j)
    {
      for(int i=0; i < n_repeat; ++i)
	{
	  Enemin[j][i] = 0;
	}
    }
  
  // temp
  ntmax[0] = 0;
  ntmax[1] = static_cast<int>(ts[3]/dt);
  ntmax[2] = static_cast<int>(ts[4]/dt);
  ntmax[3] = static_cast<int>(ts[5]/dt);
  ntmax[4] = static_cast<int>(ts[8]/dt);
  // Qs_f, Qs_i => Init()

  sprintf(dname[eDat], "dir-dat");
  sprintf(dname[eEtc], "dir-etc");
  // sinc()--------------------------------------------------
  // st
  sprintf(fname_lst[eGRAPH_ST_R0][eSinc   ][eName], "%s/%s_sinc_st.dat"   , dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_ST_R1][eSinc   ][eName], "%s/%s_r1_sinc_st.dat", dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_ST_R2][eSinc   ][eName], "%s/%s_r2_sinc_st.dat", dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_ST_R3][eSinc   ][eName], "%s/%s_r3_sinc_st.dat", dname[eEtc], name.c_str());
  // ev
  sprintf(fname_lst[eGRAPH_EV_R0][eSinc   ][eName], "%s/%s_sinc_ev.dat"   , dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_EV_R1][eSinc   ][eName], "%s/%s_r1_sinc_ev.dat", dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_EV_R2][eSinc   ][eName], "%s/%s_r2_sinc_ev.dat", dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_EV_R3][eSinc   ][eName], "%s/%s_r3_sinc_ev.dat", dname[eEtc], name.c_str());
  // optimal_result()--------------------------------------------------
  // st
  sprintf(fname_lst[eGRAPH_ST_R0][eSaiteki][eTemp], "%s/%s_saiteki_st_%%05d.dat"   , dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_ST_R1][eSaiteki][eTemp], "%s/%s_r1_saiteki_st_%%05d.dat", dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_ST_R2][eSaiteki][eTemp], "%s/%s_r2_saiteki_st_%%05d.dat", dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_ST_R3][eSaiteki][eTemp], "%s/%s_r3_saiteki_st_%%05d.dat", dname[eEtc], name.c_str());
  // ev
  sprintf(fname_lst[eGRAPH_EV_R0][eSaiteki][eTemp], "%s/%s_saiteki_ev_%%05d.dat"   , dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_EV_R1][eSaiteki][eTemp], "%s/%s_r1_saiteki_ev_%%05d.dat", dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_EV_R2][eSaiteki][eTemp], "%s/%s_r2_saiteki_ev_%%05d.dat", dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_EV_R3][eSaiteki][eTemp], "%s/%s_r3_saiteki_ev_%%05d.dat", dname[eEtc], name.c_str());
  
  for(int i = eGRAPH_ST_R1; i <= eENERGY; ++i)
    {// non use eGRAPH_ST_R0, eGRAPH_EV_R0 
      sprintf(fname_lst[i][eSaiteki][eName],
	      fname_lst[i][eSaiteki][eTemp], l_nnn );
      printf("[%2d]\n", i);
      printf("eTemp:%p:%s\n", fname_lst[i][eSaiteki][eTemp], fname_lst[i][eSaiteki][eTemp]);
      printf("eName:%p:%s\n", fname_lst[i][eSaiteki][eName], fname_lst[i][eSaiteki][eName]);
    }
  // 1.12
  
  // --------------------------------------------------
  // improvement rate
  sprintf(fname_lst[eENEMIN2][eSaiteki][eName], "%s/%s_enemin2.dat", dname[eEtc], name.c_str());
  // --------------------------------------------------
  sprintf(fname_lst[eENERGY][eSinc   ][eName], "%s/%s_energy_sinc.dat", dname[eEtc], name.c_str());
  sprintf(fname_lst[eENERGY][eSaiteki][eTemp], "%s/%s_energy_saiteki_%%05d.dat", dname[eEtc], name.c_str());

  // zu()
  sprintf(fname_lst[eZU_E1][eSinc   ][eTemp], "%s/%s_sinc_e1_%%05d.dat"   , dname[eDat], name.c_str()); // 
  sprintf(fname_lst[eZU_E2][eSinc   ][eTemp], "%s/%s_sinc_e2_%%05d.dat"   , dname[eDat], name.c_str()); // 
  sprintf(fname_lst[eZU_E3][eSinc   ][eTemp], "%s/%s_sinc_e3_%%05d.dat"   , dname[eDat], name.c_str()); // 
  sprintf(fname_lst[eZU_E1][eSaiteki][eTemp], "%s/%s_saiteki_e1_%%05d.dat", dname[eDat], name.c_str() );
  sprintf(fname_lst[eZU_E2][eSaiteki][eTemp], "%s/%s_saiteki_e2_%%05d.dat", dname[eDat], name.c_str() );
  sprintf(fname_lst[eZU_E3][eSaiteki][eTemp], "%s/%s_saiteki_e3_%%05d.dat", dname[eDat], name.c_str() );

  // trajectory
  sprintf(fname_lst[eGRAPH_ENDEFFECTOR][eSinc][eName]
	  ,"%s/%s_sinc_endeffector.dat"
	  , dname[eEtc], name.c_str());
  sprintf(fname_lst[eGRAPH_ENDEFFECTOR][eSaiteki][eTemp]
	  ,"%s/%s_saiteki_endeffector_%%05d.dat"
	  , dname[eEtc], name.c_str());
  
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


  for(uint i=0; i < eOutputFileNameSize; ++i){
    for(uint j=0; j < eWhereSize; ++j){
      printf("[%2d]\t", i);
      printf("Name=%s\t", fname_lst[i][j][eName]);
      printf("Temp=%s\n", fname_lst[i][j][eTemp]);
    }
  }

  //for(int i = eZU_E1; i < eOutputFileNameSize; ++i)
  //  {
  //    printf("[%2d]\n", i);
  //    printf("eTemp:%p:%s\n", fname_lst[i][eSinc][eTemp], fname_lst[i][eSinc][eTemp]);
  //    printf("eName:%p:%s\n", fname_lst[i][eSinc][eName], fname_lst[i][eSinc][eName]);
  //    printf("eTemp:%p:%s\n", fname_lst[i][eSaiteki][eTemp], fname_lst[i][eSaiteki][eTemp]);
  //    printf("eName:%p:%s\n", fname_lst[i][eSaiteki][eName], fname_lst[i][eSaiteki][eName]);
  //  }
  //  assert(false);

  //  assert(false);
}

IDP::~IDP()
{//  delete Enemin;
}

void
IDP::Init(Model _state_i,//ScaraRobot _state_i,VRManipulator _v1_i,
	  SinCurv _sc,
	  double _Qs_i, double _Qs_f,
	  double _Qv_i, double _Qv_f,
	  double _Qa_i, double _Qa_f )
{
  //Qs_f   = ; => initialize at Inheritance
  // Qs_f ? Sth[n_repeat]
  state_i     = _state_i;
  state_w     = _state_i;		// now working states
  keirosokudo =      0.0;
  sc	      = _sc     ;
  Qs_i = _Qs_i;  Qs_f = _Qs_f;
  Qv_i = _Qv_i;  Qv_f = _Qv_f;
  Qa_i = _Qa_i;  Qa_f = _Qa_f;
  
  for(int i=0; i <= n_repeat - 2 ; ++i)
    {
      Sth[i] = sc.sin_th1(t1*i);
    }

  for(uint i=0; i < sizeof(make_animation_flag)/sizeof(*make_animation_flag); ++i)
    {
      make_animation_flag[i] = 1;
    }
}

void
IDP::start()
{
  debug();
  sinc();			// 初期軌道を与える
  
  state_w = state_i;
  enemin2 = 0.0;
  
  // printf("enemin2=%lf\n", enemin2);
  // printf("e1.energy=%lf\n", e1.energy);
  // printf("e2.energy=%lf\n", e2.energy);
  // printf("e3.energy=%lf\n", e3.energy);
  // assert(false);
  // ofstream ofs(fname_enemin2);
  ofstream ofs(fname_lst[eENEMIN2][eSaiteki][eName]);
  for(l_nnn = 1; l_nnn <= n_repeat; ++l_nnn)
    {
      repeat();
      keirosokudo = 0.0;
      tansaku();
      //Enemin[l_nnn] = optimal_result(0);
      //Enemin[l_nnn] = optimal_result( make_animation_flag[l_nnn] );
      optimal_result( make_animation_flag[l_nnn] );
      Enemin[0][l_nnn] = e1.energy + e2.energy + e3.energy;
      Enemin[1][l_nnn] = e1.energy;
      Enemin[2][l_nnn] = e2.energy;
      Enemin[3][l_nnn] = e3.energy;
      
      ofs << l_nnn            <<" "
	  << Enemin[0][l_nnn] <<" " // [r1 + r2 + r3]
	  << Enemin[1][l_nnn] <<" " // [r1]
	  << Enemin[2][l_nnn] <<" " // [r2]
	  << Enemin[3][l_nnn] <<" " // [r3]
	  << kiza
	  << endl;
      
      cout << name  <<" "
	   << l_nnn <<" "
	   << Enemin[0][l_nnn] <<" " // [r1 + r2 + r3]
	   << Enemin[1][l_nnn] <<" " // [r1]
	   << Enemin[2][l_nnn] <<" " // [r2]
	   << Enemin[3][l_nnn] <<" " // [r3]
	   << kiza << endl;
      
      if ( l_nnn > 15 && Enemin[0][l_nnn] > Enemin[0][l_nnn - 1] )
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
  for(nt = 1; nt <= ntmax; nt++)
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
  // rewrite Functa 
  //cout << __FILE__ <<"-"<<__FUNCTION__<<"-"<< __LINE__ << endl;
  set_position();
  GraphicalObject *obj[] =
    {
      new Circle3D( r1.x0, r1.y0, r1.z0, 0.01),
      new Circle3D( r1.x , r1.y , r1.z0, 0.01),
      new Circle3D( r2.x0, r2.y0, r2.z0, 0.01),
      new Circle3D( r2.x , r2.y , r2.z0, 0.01),
      new Circle3D( r3.x0, r3.y0, r3.z0, 0.01),
      new Circle3D( r3.x0, r3.y0, r3.z , 0.01),
      new Line3D(   r1.x0, r1.y0, r1.z0, r1.x , r1.y , r1.z0 ),
      new Line3D(   r2.x0, r2.y0, r2.z0, r2.x , r2.y , r2.z0 ),
      new Line3D(   r3.x0, r3.y0, r3.z0, r3.x0, r3.y0, r3.z  ),
      new Circle3D( v1.x0, v1.y0, v1.z0, 0.01),
      new Circle3D( v1.x , v1.y , v1.z0, 0.01),
      new Line3D(   v1.x0, v1.y0, v1.z0, v1.x , v1.y , v1.z0 ),
    };
  
  //  Circle3D real2_p2(r3.x0, r3.y0, r3.z0, 0.01);
  //  Line3D real2_ln(r2.x0, r2.y0, r1.z0, r2.x,  r2.y , r1.z0 );
  
  for(uint i = 0; i < sizeof(obj)/sizeof(*obj); ++i)
    {
      obj[i]->draw(ofs);
      delete obj[i];
    }
  // Functa
  // function pointer
  // graph_*(fp)
}

void IDP::zu(ostream& ofs1, ostream& ofs2)
{
  //cout << __FILE__ <<"-"<<__FUNCTION__<<"-"<< __LINE__ << endl;
  set_position();
  
  GraphicalObject *obj[] =
    {
      // ofs1
      new Circle3D(r1.x0, r1.y0, r1.z0, 0.01), // 1
      new Circle3D(r1.x , r1.y , r1.z0, 0.01), // 2
      new Circle3D(r2.x0, r2.y0, r2.z0, 0.01), // 3
      new Circle3D(r2.x , r2.y , r2.z0, 0.01), // 4
      new Circle3D(r3.x0, r3.y0, r3.z0, 0.01), // 5
      new Circle3D(r3.x0, r3.y0, r3.z , 0.01), // 6
      new Line3D(  r1.x0, r1.y0, r1.z0, r1.x,  r1.y , r1.z0 ), // 7
      new Line3D(  r2.x0, r2.y0, r2.z0, r2.x,  r2.y , r2.z0 ), // 8
      new Line3D(  r3.x0, r3.y0, r3.z0, r3.x0, r3.y0, r3.z  ), // 9
      // ofs2
      new Circle3D(v1.x0, v1.y0, v1.z0, 0.01), // 10
      new Circle3D(v1.x , v1.y , v1.z0, 0.01), // 11
      new Line3D(  v1.x0, v1.y0, v1.z0,v1.x , v1.y , v1.z0 ), // 12
    };

  for(uint i = 0; i < sizeof(obj)/sizeof(*obj); ++i)
    {
      if (i < 9)
	obj[i]->draw(ofs1);
      else
	obj[i]->draw(ofs2);
      delete obj[i];
    }
}

void IDP::zu(ostream& ofs1, ostream& ofs2, ostream& ofs3)
{
  //cout << __FILE__ <<"-"<<__FUNCTION__<<"-"<< __LINE__ << endl;
  set_position();
  GraphicalObject *obj[] =
    {// ofs1
      new Circle3D(r1.x0, r1.y0, r1.z0, 0.01), // 1
      new Circle3D(r1.x , r1.y , r1.z0, 0.01), // 2
      new Circle3D(r2.x0, r2.y0, r2.z0, 0.01), // 3
      new Circle3D(r2.x , r2.y , r2.z0, 0.01), // 4
      new Circle3D(r3.x0, r3.y0, r3.z0, 0.01), // 5
      new Circle3D(r3.x0, r3.y0, r3.z , 0.01), // 6
      new Line3D(  r1.x0, r1.y0, r1.z0, r1.x , r1.y , r1.z0 ), // 7
      new Line3D(  r2.x0, r2.y0, r2.z0, r2.x , r2.y , r2.z0 ), // 8
      new Line3D(  r3.x0, r3.y0, r3.z0, r3.x0, r3.y0, r3.z  ), // 9
      // ofs2
      new Circle3D(v1.x0, v1.y0, v1.z0, 0.01), // 10
      new Circle3D(v1.x , v1.y , v1.z0, 0.01), // 11
      new Line3D(v1.x0, v1.y0, v1.z0,v1.x , v1.y , v1.z0 ), // 12
      // ofs3
      new Line3D(pnt_i.x, pnt_i.y, pnt_i.z, pnt_c.x, pnt_c.y, pnt_c.z),
    };

  for(uint i = 0; i < sizeof(obj)/sizeof(*obj); ++i)
    {
      if (i < 9)
	obj[i]->draw(ofs1);
      else if (i >= 9 && i < 12)
	obj[i]->draw(ofs2);
      else
	obj[i]->draw(ofs3);
      delete obj[i];
    }
}

void IDP::zu(ostream& ofs1, ostream& ofs2, ostream& ofs3, vector<Point> trjc)
{
  //cout << __FILE__ <<"-"<<__FUNCTION__<<"-"<< __LINE__ << endl;
  set_position();
  GraphicalObject *obj[] =
    {// ofs1
      new Circle3D(r1.x0, r1.y0, r1.z0, 0.01), // 1
      new Circle3D(r1.x , r1.y , r1.z0, 0.01), // 2
      new Circle3D(r2.x0, r2.y0, r2.z0, 0.01), // 3
      new Circle3D(r2.x , r2.y , r2.z0, 0.01), // 4
      new Circle3D(r3.x0, r3.y0, r3.z0, 0.01), // 5
      new Circle3D(r3.x0, r3.y0, r3.z , 0.01), // 6
      new Line3D(  r1.x0, r1.y0, r1.z0, r1.x,  r1.y , r1.z0 ), // 7
      new Line3D(  r2.x0, r2.y0, r2.z0, r2.x,  r2.y , r2.z0 ), // 8
      new Line3D(  r3.x0, r3.y0, r3.z0, r3.x0, r3.y0, r3.z  ), // 9
      // ofs2
      new Circle3D(v1.x0, v1.y0, v1.z0, 0.01), // 10
      new Circle3D(v1.x , v1.y , v1.z0, 0.01), // 11
      new Line3D(  v1.x0, v1.y0, v1.z0, v1.x , v1.y , v1.z0 ), // 12
      // ofs3
      //new Line3D(pnt_i.x, pnt_i.y, pnt_i.z, pnt_c.x, pnt_c.y, pnt_c.z),
      //new Line3D(pnt_i.x, pnt_i.y, pnt_i.z, pnt_c.x, pnt_c.y, pnt_c.z),
      //new Line3D(pnt_c.x, pnt_c.y, pnt_i.z, pnt_f.x, pnt_f.y, pnt_c.z),
    };
  
  for(uint i = 0; i < sizeof(obj)/sizeof(*obj); ++i)
    {
      if (i < 9)
	obj[i]->draw(ofs1);
      else if (i >= 9 && i < 12)
	obj[i]->draw(ofs2);
    }

  for(uint i = 1; i < trjc.size(); ++i)
    {
      Line3D line(trjc[i-1].x, trjc[i-1].y, trjc[i-1].z,
		  trjc[i  ].x, trjc[i  ].y, trjc[i  ].z);
      line.draw(ofs3);
    }
}

void IDP::makeGifFile() const
{
  vector<string> fnames;
  char fname_temp[64];
  fnames.push_back(fname_lst[eENERGY][eSinc][eName]);  
  for(uint i = 1;
      i < sizeof(make_animation_flag)/sizeof(*make_animation_flag);
      ++i)
    {
      if ( make_animation_flag[i] )
	{
	  //char fname[32]; sprintf(fname, ftemp_energy, i);
	  //fnames.push_back(fname);
	  sprintf( fname_temp,
		   fname_lst[eENERGY][eSaiteki][eTemp], i );
	  fnames.push_back( fname_temp );
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
  fprintf(gp, ",'%s' u 1:2 w l ls 3 title 'result'\n"
	  , fnames.back().c_str() );  

  pclose(gp);
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
	    at[i]   = Mat[i][at[i+1] ][at[i+2] ];
	    Qs[i]   = Sth[i] + kiza*(at[i] - N3);
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
    at[ii]   = Mat[ii][at[ii+1] ][at[ii+2] ];
    sat[ii]  =  at[ii];
  }  
  //  for(int ii = 6; ii > 0; ii--){//Magick
  for(int ii = n_divide-2; ii > 0; --ii){//Magick
    Qs[ii] = Sth[ii] + kiza*(at[ii] - N3);
  }
  
  //cur = 4; ((8 - 2) - 3) + 1//n_divide=8
  cur = ((n_divide - 2) - 3) + 1;
  settei(cur);
  //    reinitialize();
  
  sprintf(fname_lst[eENERGY][eSaiteki][eName],
	  fname_lst[eENERGY][eSaiteki][eTemp], l_nnn);
  
  sprintf(fname_lst[eGRAPH_ENDEFFECTOR][eSaiteki][eName],
	  fname_lst[eGRAPH_ENDEFFECTOR][eSaiteki][eTemp], l_nnn);
  
  ofstream ofs_st   , ofs_ev   ;
  ofstream ofs_r1_st, ofs_r1_ev;
  ofstream ofs_r2_st, ofs_r2_ev;
  ofstream ofs_r3_st, ofs_r3_ev;
  ofstream ofs_energy;
  ofstream ofs_endeffector;

  FILE *fp_lst[eOutputFileNameSize];

  if (make_animation_flag)
    {
      printf("%s-%d\n", __FUNCTION__, __LINE__);
      for(int i = eGRAPH_ST_R1; i <= eENERGY; ++i)
	{// non use eGRAPH_ST_R0, eGRAPH_EV_R0 
	  sprintf(fname_lst[i][eSaiteki][eName],
		  fname_lst[i][eSaiteki][eTemp], l_nnn       );
	  printf("eName:%s\t", fname_lst[i][eSaiteki][eName] );
	  printf("eTemp:%s\n", fname_lst[i][eSaiteki][eTemp] );
	  fp_lst[i] = fopen(fname_lst[i][eSaiteki][eName],"w");
	}
    }
  //  assert(false);
  
  //  for(int nt = 1; nt <= ntmax[ cur ]; nt++)
  vector<Point>trjc;
  
  state_w = state_i;
  enemin2 = 0.0;
  for(int nt = 1; nt <= ntmax[ cur ]; nt++)
    {
      t += dt;
      sim(cur, t);
      enemin2 = En();

      if (make_animation_flag)
	{
	  graph_standard(fp_lst[eGRAPH_ST_R1],
			 t    ,
			 r1.th, r1.thv, r1.tha,
			 e1.ev, e1.ia , e1.tau,
			 e1.energy    , e1.ev*e1.ia);
	  
	  graph_volt_constituent(fp_lst[eGRAPH_EV_R1],
				 t, e1.ev ,
				 b1*r1.thv, b2*r1.tha, b3*e1.tau);
	  // --------------------------------------------------
	  graph_standard(fp_lst[eGRAPH_ST_R2],
			 t    ,
			 r2.th, r2.thv, r2.tha,
			 e2.ev, e2.ia , e2.tau,
			 e2.energy    , e2.ev*e2.ia);
	  graph_volt_constituent(fp_lst[eGRAPH_EV_R2],
				 t, e2.ev ,
				 b1*r2.thv, b2*r2.tha, b3*e2.tau);
	  // --------------------------------------------------
	  graph_standard(fp_lst[eGRAPH_ST_R3],
			 t    ,
			 // r3.z , r3.zv, r3.za ,
			 r3.z0 , r3.zv, r3.za ,			 
			 e3.ev, e3.ia, e3.tau,
			 e3.energy, e3.ev*e3.ia);
	  graph_volt_constituent(fp_lst[eGRAPH_EV_R3],
				 t    ,
				 e3.ev,
				 b1*r3.Thv(), b2*r3.Tha(), b3*e3.tau);
	  // --------------------------------------------------	  
	  //graph_endeffector(ofs_endeffector,  t, v1, r3);
	  //fprintf(fp_foo[eEnergy], "%lf %lf %lf %lf %lf\n"
	  graph_endeffector(fp_lst[eGRAPH_ENDEFFECTOR], t, v1, r3);
	  fprintf(fp_lst[eENERGY],
		  "%lf %lf %lf %lf %lf\n"
		  , t
		  , enemin2	// [r1+r2+r3]
		  , e1.energy	// [r1]
		  , e2.energy	// [r2]
		  , e3.energy ); //[r3]
	}
      
      if (nt % 10 == 0)
	{
	  trjc.push_back( Point(v1.x, v1.y, r3.z0) );
	  // enum => define
	  // like variable.h
	  sprintf(fname_lst[eZU_E1][eSaiteki][eName],
		  fname_lst[eZU_E1][eSaiteki][eTemp], nt);
	  sprintf(fname_lst[eZU_E2][eSaiteki][eName],
		  fname_lst[eZU_E2][eSaiteki][eTemp], nt);
	  sprintf(fname_lst[eZU_E3][eSaiteki][eName],
		  fname_lst[eZU_E3][eSaiteki][eTemp], nt);
	  
	  ofstream ofs1(fname_lst[eZU_E1][eSaiteki][eName]);
	  ofstream ofs2(fname_lst[eZU_E2][eSaiteki][eName]);
	  ofstream ofs3(fname_lst[eZU_E3][eSaiteki][eName]);
	  
	  zu(ofs1, ofs2, ofs3, trjc);
	}
    } //step nt
  
  if ( make_animation_flag )
    {
      for(int i = eGRAPH_ST_R1; i < eENERGY; ++i)
	{
	  fclose( fp_lst[i] );
	}
    }
  
  return enemin2;
}

// IDP::settei => inline IDP::settei()
// 
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

// inline IDP::sim()
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
      //r3.z0 = th1[2] + th1v[2]*(t - ts[2]) + 0.5*r3.za*pow((t - ts[2]), 2);
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
      //r3.z0  = th1[4] + th1v[4]*(t - ts[4]) + 0.5*r3.za*pow((t - ts[4]), 2);
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
      //r3.z0 = th1[6] + th1v[6]*(t - ts[6]) + 0.5*r3.za*pow((t - ts[6]), 2);
    }
  if (ts[7] < t && t < ts[8])
    {
      r3.za = Qa[8];
      r3.zv = Qv[7] + r3.za*(t - ts[7]);
      r3.z0 = Qs[7] + Qv[7]*(t - ts[7]) + 0.5*r3.za*SQ(t - ts[7]);
      //r3.z0  = th1[7] + th1v[7]*(t - ts[7]) + 0.5*r3.za*pow((t - ts[7]), 2);
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

void IDP::set_position(){}

Model::Model(){}
Model::Model(VRManipulator _v1,
	     Revolute _r1     , Revolute _r2, Prismatic _r3,
	     Energy   _e1     , Energy   _e2, Energy    _e3 )
  : v1(_v1),
    r1(_r1), r2(_r2), r3(_r3),
    e1(_e1), e2(_e2), e3(_e3){}
Model::~Model(){}

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
#include <iterator>
#include "Manipulator.hpp"
#include "IDP.hpp"
#include "Bi3_2.hpp"

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
	 double _t_i, double _t_f, double _dt,
	 double _kiza_x, double _kiza_y, double _kiza_z,
	 int _n_repeat, int _n_divide, int _NN)
  : name(_name), t_i(_t_i), t_f(_t_f), dt(_dt),
    kiza_x(_kiza_x), kiza_y(_kiza_y), kiza_z(_kiza_z),
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

  //  INITARY(Sth);
  INITARY(Sxx);
  INITARY(Syy);
  INITARY(Szz);
  INITARY(ntmax);
  //INITARY(Qs); //INITARY(Qv); //INITARY(Qa);
  //INITARY(at); //INITARY(sat);
  INITARY(Xs); INITARY(Xv); INITARY(Xa);
  INITARY(Ys); INITARY(Yv); INITARY(Ya);
  INITARY(Zs); INITARY(Zv); INITARY(Za);
  
  for(int i = 0; i <= n_divide - 2; ++i)
    {//sat[i] = N3;
      sxt[i] = N3;
      syt[i] = N3;
      szt[i] = N3;
    }

  for(int cur = 0; cur <= n_divide - 2; ++cur)
    for(int nx1 = 0; nx1 <= NN; ++nx1)
      for(int ny1 = 0; ny1 <= NN; ++ny1)
	for(int nx2 = 0; nx2 <= NN; ++nx2)
	  for(int ny2 = 0; ny2 <= NN; ++ny2)
	    {// Mat[i][j][k] = 0;
	      Mxt[cur][nx1][ny1][nx2][ny2] = 0;
	      Myt[cur][nx1][ny1][nx2][ny2] = 0;
	      //Mzt[cur][nx1][ny1][nx2][ny2] = 0;	  
	    }

  for(int at1 = 1; at1 <= NN; at1++)
    for(int at2 = 1; at2 <= NN; at2++)
      for(int at3 = 1; at3 <= NN; at3++)
	{
	  Energ[at1][at2][at3] = 0.0;
	}
  
  for(int nx1=0; nx1 <= NN; ++nx1)
    for(int ny1=0; ny1 <= NN; ++ny1)
      for(int nx2=0; nx2 <= NN; ++nx2)
	for(int ny2=0; ny2 <= NN; ++ny2)
	  for(int nx3=0; nx3 <= NN; ++nx3)
	    for(int ny3=0; ny3 <= NN; ++ny3)
	      {		// n4, n5, n6
		Energ2[nx1][ny1][nx2][ny2][nx3][ny3] = 0.0;
	      }

  for(int nx1=0; nx1 <= NN; ++nx1)
    for(int ny1=0; ny1 <= NN; ++ny1)
      for(int nz1=0; nz1 <= NN; ++nz1)
	for(int nx2=0; nx2 <= NN; ++nx2)
	  for(int ny2=0; ny2 <= NN; ++ny2)
	    for(int nz2=0; nz2 <= NN; ++nz2)
	      for(int nx3=0; nx3 <= NN; ++nx3)
		for(int ny3=0; ny3 <= NN; ++ny3)
		  for(int nz3=0; nz3 <= NN; ++nz3)
		    {		// n4, n5, n6
		      //Energ3[nx1][ny1][nz1][nx2][ny2][nz2][nx3][ny3][nz3];
		    }
    
  //Enemin = new double(n_repeat); // hang in the air. offset
  for(int i=0; i < n_repeat; ++i)
    {
      Enemin[i] = 0;
    }
  
  // temp
  ntmax[0] = 0;
  ntmax[1] = static_cast<int>(ts[3]/dt);
  ntmax[2] = static_cast<int>(ts[4]/dt);
  ntmax[3] = static_cast<int>(ts[5]/dt);
  ntmax[4] = static_cast<int>(ts[8]/dt);
  // Qs_f, Qs_i => Init()

  // --------------------------------------------------
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

  // energy
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

  // r1, r2, sinc_st, sinc_ev
  sprintf(fname_r1_sinc_st, "%s/%s_r1_sinc_st.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r2_sinc_st, "%s/%s_r2_sinc_st.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r3_sinc_st, "%s/%s_r3_sinc_st.dat"
	  , dname_dir_etc, name.c_str() );
  
  sprintf(fname_r1_sinc_ev, "%s/%s_r1_sinc_ev.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r2_sinc_ev, "%s/%s_r2_sinc_ev.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r3_sinc_ev, "%s/%s_r3_sinc_ev.dat"
	  , dname_dir_etc, name.c_str() );

  // r1, r2 saiteki_st, sinc_ev
  sprintf(ftemp_r1_saiteki_st, "%s/%s_r1_saiteki_st_%%05d.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(ftemp_r2_saiteki_st, "%s/%s_r2_saiteki_st_%%05d.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(ftemp_r3_saiteki_st, "%s/%s_r3_saiteki_st_%%05d.dat"
	  , dname_dir_etc, name.c_str() );

  sprintf(ftemp_r1_saiteki_ev, "%s/%s_r1_saiteki_ev_%%05d.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(ftemp_r2_saiteki_ev, "%s/%s_r2_saiteki_ev_%%05d.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(ftemp_r3_saiteki_ev, "%s/%s_r3_saiteki_ev_%%05d.dat"
	  , dname_dir_etc, name.c_str() );

  sprintf(fname_r1_enemin2, "%s/%s_r1_enemin2.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r2_enemin2, "%s/%s_r2_enemin2.dat"
	  , dname_dir_etc, name.c_str() );
  sprintf(fname_r3_enemin2, "%s/%s_r3_enemin2.dat"
	  , dname_dir_etc, name.c_str() );


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
}

IDP::~IDP()
{
  //  delete Enemin;
}

// void
// IDP::Init(Model	_state_i,	//ScaraRobot _state_i,VRManipulator _v1_i,
// 	  SinCurv _sc,
// 	  double _Qs_i, double _Qs_f,
// 	  double _Qv_i, double _Qv_f,
// 	  double _Qa_i, double _Qa_f)
// {
//   state_i  = _state_i;
//   state_w  = _state_i;		// now working states
  
//   sc	   = _sc;
//   Qs_i = _Qs_i;  Qs_f = _Qs_f;
//   Qv_i = _Qv_i;  Qv_f = _Qv_f;
//   Qa_i = _Qa_i;  Qa_f = _Qa_f;
  
//   //  for(int i=0; i <= n_repeat - 2 ; ++i)
//   for(int i=0; i <= n_divide - 2 ; ++i)  
//     {
//       Sth[i] = sc.sin_th1(t1*i);
//     }

//}

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
      
      ofs << l_nnn <<" "
	  << Enemin[l_nnn] <<" "
	  << kiza_x <<" "
	  << kiza_y <<" "
	  << kiza_z <<" "
	  << endl;
      
      //cout << name <<" "<<l_nnn <<" "<<Enemin[l_nnn]<<" "<< kiza << endl;
      if ((l_nnn > 15) && (Enemin[l_nnn] > Enemin[l_nnn - 1]))
	{
	  //kiza *= 1.2;
	  kiza_x *= 1.2; kiza_y *= 1.2; kiza_z *= 1.2;
	}

      cout << "l_nnn = " << l_nnn << endl;

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
  //Energ[at[cur] ][at[cur+1] ][at[cur+2] ] = energy;//=> sum_ev_ia_dt;
  Energ2[xt[cur]][yt[cur]][xt[cur+1]][yt[cur+1]][xt[cur+2]][yt[cur+2]]
    = energy;//=> sum_ev_ia_dt;
  //Energ3[at[cur] ][at[cur+1] ][at[cur+2] ] = energy;//=> sum_ev_ia_dt;
}

void IDP::makeGifFile() const
{
  vector<string> fnames;
  for(uint i = 0;
      i < sizeof(make_animation_flag)/sizeof(*make_animation_flag);
      ++i)
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
  fprintf(gp, "'%s' u 1:2 w l ls 1 title 'sinc'"
	  , fnames[0].c_str() );
  for(uint i = 1; i < fnames.size() - 1; ++i)
    {
      fprintf(gp, ",'%s' u 1:2 w l ls 2 title '%d' "
	      , fnames[i].c_str(), i);
    }
  fprintf(gp, ",'%s' u 1:2 w l ls 3 title 'result'\n"
	  , fnames.back().c_str() );  

  pclose(gp);
}

void IDP::repeat()		// ofstream
{
  //int ff1;
  int ffx, ffy, ffz;
  //Magick
  for(int ii=6; ii > 0; --ii)
    {  /* 中心位置に探索結果の最適経路を選ぶ */
      //Sth[ii] = Sth[ii] + kiza*(sat[ii] - N3);
      Sxx[ii] = Sxx[ii] + kiza_x*(sxt[ii] - N3);
      Syy[ii] = Syy[ii] + kiza_y*(syt[ii] - N3);
      Szz[ii] = Szz[ii] + kiza_z*(szt[ii] - N3);
    }
  
  //{// get summation
  //  double sum = 0.0;
  //  for(int ii=6; ii > 0; --ii){//Magick
  //    sum = sum + fabs(sat[ii] - N3);
  //  }
  //  ff1 = static_cast<int>(sum);
  //}

  {// get summation
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    for(int ii=6; ii > 0; --ii){//Magick
      sum_x = sum_x + fabs(sxt[ii] - N3);
      sum_y = sum_y + fabs(syt[ii] - N3);
      sum_z = sum_z + fabs(szt[ii] - N3);
    }
    ffx = static_cast<int>(sum_x);
    ffy = static_cast<int>(sum_y);
    ffz = static_cast<int>(sum_z);
  }
  
  if (0 == ffx)
    {
      kiza_x = kiza_x*0.8;
      printf("%s-%s, kiza_x=%lf, nnn=%d\n"
	     , __FILE__, __FUNCTION__, kiza_x, l_nnn);
    }
  if (0 == ffy)
    {
      kiza_y = kiza_y*0.8;
      printf("%s-%s, kiza_y=%lf, nnn=%d\n"
	     , __FILE__, __FUNCTION__, kiza_y, l_nnn);
    }
  if (0 == ffz)
    {
      kiza_z = kiza_z*0.8;
      printf("%s-%s, kiza_z=%lf, nnn=%d\n"
	     , __FILE__, __FUNCTION__, kiza_z, l_nnn);
    }
}

void IDP::tansaku()
{
  // at[cur] = 1 : NN;
  // cur = 1 : (n_divide - 2) - 2;
  for(int cur = 1; cur <= 4; cur++){//Magick latter, replace the divide number
    //--------------------------------------------------
    for(xt[cur]=1; xt[cur] <= NN; xt[cur]++){
      for(yt[cur]=1; yt[cur] <= NN; yt[cur]++){
	//--------------------------------------------------
	for(xt[cur+1]=1; xt[cur+1] <= NN; xt[cur+1]++){
	  for(yt[cur+1] = 1; yt[cur+1] <= NN; yt[cur+1]++){
	    // --------------------------------------------------
	    for(xt[cur+2] = 1; xt[cur+2] <= NN; xt[cur+2]++){
	      for(yt[cur+2]=1; yt[cur+2] <= NN; yt[cur+2]++){
		// --------------------------------------------------
		//Qs[cur+2] = Sth[cur+2] + kiza*(at[cur+2] - N3);	  
		//Qs[cur+1] = Sth[cur+1] + kiza*(at[cur+1] - N3);
		//Qs[cur  ] = Sth[cur  ] + kiza*(at[cur  ] - N3);
		Xs[cur+2] = Sxx[cur+2] + kiza_x*(xt[cur+2] - N3);
		Xs[cur+1] = Sxx[cur+1] + kiza_x*(xt[cur+1] - N3);
		Xs[cur+0] = Sxx[cur+0] + kiza_x*(xt[cur+0] - N3);
		Ys[cur+2] = Syy[cur+2] + kiza_y*(yt[cur+2] - N3);
		Ys[cur+1] = Syy[cur+1] + kiza_y*(yt[cur+1] - N3);
		Ys[cur+0] = Syy[cur+0] + kiza_y*(yt[cur+0] - N3);
		// Zs[cur+2] = Szz[cur+2] + kiza_z*(szt[cur+2] - N3);
		// Zs[cur+1] = Szz[cur+1] + kiza_z*(szt[cur+1] - N3);
		// Zs[cur+0] = Szz[cur+0] + kiza_z*(szt[cur+0] - N3);
		for(int i = cur; --i > 0;)
		  {
		    //at[i]  = Mat[i][at[i+1] ][at[i+2] ];
		    //Qs[i] = Sth[i] + kiza*(at[i] - N3);
		    xt[i] = Mxt[i][xt[i+1]][yt[i+1]][xt[i+2]][yt[i+2]];
		    yt[i] = Myt[i][xt[i+1]][yt[i+1]][xt[i+2]][yt[i+2]];
		    //zt[i] = Mzt[i][zt[i+1] ][zt[i+2] ];
		    Xs[i] = Sxx[i] + kiza_x*(xt[i] - N3);
		    Ys[i] = Syy[i] + kiza_y*(yt[i] - N3);
		    //Zs[i] = Szz[i] + kiza_z*(zt[i] - N3);
		  }
		settei(cur);// set  [XYZ]v\[\], [XYZ]a\[\];
		En_sim( ntmax[cur], cur );
		// 
	      }
	    }
	    // [cur + 0]
	  }
	}
	// [cur + 1]
      }
    } // [cur + 2]
    //Hikaku(cur);
    Hikaku2(cur);
    //Hikaku3(cur);
  }
  return;
}

double
IDP::optimal_result(int make_animation_flag)
{
  int cur;
  double t = 0;

  //at[6] = sat[6];
  //at[5] = sat[5];
  //at[4] = sat[4];
  for(int i = 6; i > 3; i--)
    {
      xt[i] = sxt[i];
      yt[i] = syt[i];
      //zt[i] = szt[i];
    }
  //  for(int ii = ((n_divide-2)-3); ii > 0; --ii){//Magick
  for(int i = 3; i > 0; i--)
    {
      //at[i] = Mat[i][at[i+1] ][at[i+2] ];
      //sat[i] =  at[i];
      xt[i] = Mxt[i][ xt[i+1] ][ yt[i+1] ][ xt[i+2] ][ yt[i+2] ];
      yt[i] = Myt[i][ xt[i+1] ][ yt[i+1] ][ xt[i+2] ][ yt[i+2] ];
      //zt[i] = Mzt[i][ zt[i+1] ][ zt[i+2] ];
      sxt[i] = xt[i];
      syt[i] = yt[i];
      szt[i] = zt[i];
    }
  
  //  for(int i = n_divide-2; i > 0; --i){//Magick
  for(int i = 6; i > 0; i--)
    {//Magick
      //Qs[ii] = Sth[ii] + kiza*(at[ii] - N3);
      Xs[i] = Sxx[i] + kiza_x*(xt[i] - N3);
      Ys[i] = Syy[i] + kiza_y*(yt[i] - N3);
      //Zs[i] = Szz[i] + kiza_z*(zt[i] - N3);      
    }
  
  cur = 4;// ((8 - 2) - 3) + 1//n_divide=8
  settei(cur);//cur = ((n_divide - 2) - 3) + 1;
  //    reinitialize();
  ofstream ofs_st1, ofs_ev1;//, ofs_energy1;
  ofstream ofs_st2, ofs_ev2;//, ofs_energy2;
  ofstream ofs_st3, ofs_ev3;//, ofs_energy3;

  {
    char fname[64]; sprintf(fname, ftemp_energy, l_nnn);
  }
  ofstream ofs_energy;
  if (make_animation_flag)
    {
      char fname[64];
      sprintf(fname, ftemp_r1_saiteki_st, l_nnn); ofs_st1.open(fname);
      sprintf(fname, ftemp_r1_saiteki_ev, l_nnn); ofs_ev1.open(fname);
      //sprintf(fname, ftemp_r1_energy, l_nnn); ofs_energy1.open(fname);
      
      sprintf(fname, ftemp_r2_saiteki_st, l_nnn); ofs_st2.open(fname);
      sprintf(fname, ftemp_r2_saiteki_ev, l_nnn); ofs_ev2.open(fname);
      //sprintf(fname, ftemp_energy, l_nnn); ofs_energy2.open(fname);

      sprintf(fname, ftemp_r3_saiteki_st, l_nnn); ofs_st3.open(fname);
      sprintf(fname, ftemp_r3_saiteki_ev, l_nnn); ofs_ev3.open(fname);
      //sprintf(fname, ftemp_energy, l_nnn); ofs_energy3.open(fname);
    }
  
  for(int nt = 1; nt <= ntmax[ cur ]; nt++)
    {  
      t += dt;
      sim(cur, t);
      enemin2 = En();
      
      ofs_energy << t <<" "<< enemin2 << endl;
      
      if (make_animation_flag)
	{
	  // --------------------------------------------------
	  graph_standard(ofs_st1, t,
			 r1.th, r1.thv, r1.tha,
			 e1.ev, e1.ia, e1.tau,
			 e1.energy, e1.ev*e1.ia);
	  graph_volt_constituent(ofs_ev1, t, e1.ev,
				 b1*r1.thv,
				 b2*r1.tha,
				 b3*e1.tau);
	  // --------------------------------------------------
	  graph_standard(ofs_st2, t,
			 r2.th, r2.thv, r2.tha,
			 e2.ev, e2.ia, e2.tau,
			 e2.energy, e2.ev*e2.ia);
	  graph_volt_constituent(ofs_ev2, t, e2.ev,
				 b1*r2.thv,
				 b2*r2.tha,
				 b3*e2.tau);
	  // --------------------------------------------------
	  graph_standard(ofs_st3, t,
			 r3.z, r3.zv, r3.za,
			 e3.ev, e3.ia, e3.tau,
			 e3.energy, e3.ev*e3.ia);
	  graph_volt_constituent(ofs_ev3, t, e3.ev,
				 b1*r3.Thv(),
				 b2*r3.Tha(),
				 b3*e3.tau);
	  // --------------------------------------------------
	} // fi make_animation_flag

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
	} // fi nt % 10 == 0
    } //endfor

  if (make_animation_flag)
    {
      ofs_st1.close();
      ofs_st2.close();
      ofs_st3.close();

      ofs_ev1.close();
      ofs_ev2.close();
      ofs_ev3.close();

      //ofs_energy1.close();      
      //ofs_energy2.close();
      //ofs_energy3.close();
    } // make_animation_flag
  
  return enemin2;
}

//void IDP::settei(int cur)
//{
//  //  assert(cur <= (n_divide-2)-2);
//  Qv[1] = 2.0*(Qs[1] - Qs_i)/t1;
//  Qa[1] = Qv[1]/t1;
//  Qv[2] = 2.0*(Qs[2] - Qs[1])/t1 - Qv[1];
//  Qa[2] = (Qv[2] - Qv[1])/t1;
//  Qv[3] = 2.0*(Qs[3] - Qs[2])/t1 - Qv[2];
//  Qa[3] = (Qv[3] - Qv[2])/t1;
//  if (cur == 1) return;
//    
//  Qv[4] = 2.0*(Qs[4] - Qs[3])/t1 - Qv[3];
//  Qa[4] = (Qv[4] - Qv[3])/t1;
//  if (cur == 2) return;
//    
//  Qv[5] = 2.0*(Qs[5] - Qs[4])/t1 - Qv[4];
//  Qa[5] = (Qv[5] - Qv[4])/t1;
//  if (cur == 3) return;
//    
//  Qv[6] = 2.0*(Qs[6] - Qs[5])/t1 - Qv[5];
//  Qa[6] = (Qv[6] - Qv[5])/t1;
//  // will become if (step == 4) return;
//  
//  Qs[8] = Qs_f;
//  Qv[8] = 0.0; // vary with condition
//  Qv[7] = (Qs[8] - Qs[6])/t1 - (Qv[6]+Qv[8])/2.0;
//  Qa[7] = (Qv[7] - Qv[6])/t1;
//  Qs[7] = Qs[6] + Qv[6] * t1 + Qa[7]*t1*t1/2.0;
//  Qa[8] = (Qv[8] - Qv[7])/t1;
//  if (cur == 4) return;
//
//  // if (cur == 5?)
//  //  return; ((n_divide-2)-3) + 1 < mat_max
//  // if (cur <= mat_max)
//  //     settei1()
//  // else
//  //     settei2(get the last two step acceration and velocity)
//  // 
//  // mat_max = ((n_divide - 2) - 3) + 1 (n_divide=8; 4);
//}

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

//void IDP::sim(int cur, double t)//Magick_numbers
//{
//  //  assert(cur <= (n_divide-2)-2);
//#define SQ(x) ((x)*(x))
//  if (t < ts[1])
//    {
//      r3.za = Qa[1];
//      r3.zv = r3.za*t;
//      r3.z0 = 0.5*r3.za*t*t + Qs_i;
//    }
//  if (ts[1] < t && t < ts[2])
//    {
//      r3.za = Qa[2];
//      r3.zv = Qv[1] + r3.za*(t - ts[1]);
//      r3.z0 = Qs[1] + Qv[1]*(t - ts[1]) + 0.5*r3.za*SQ(t - ts[1]);
//      //r3.z0 = th1[1] + th1v[1]*(t - ts[1]) + 0.5*r3.za*pow((t - ts[1]), 2);	
//    }
//  if (ts[2] < t && t < ts[3])
//    {
//      r3.za = Qa[3];
//      r3.zv = Qv[2] + r3.za*(t - ts[2]);
//      //	r3.z0 = th1[2] + th1v[2]*(t - ts[2]) + 0.5*r3.za*pow((t - ts[2]), 2);
//      r3.z0 = Qs[2] + Qv[2]*(t - ts[2]) + 0.5*r3.za*SQ(t - ts[2]);	
//    }
//  if (1 == cur)
//    return;
//
//  if (ts[3] < t && t < ts[4])
//    {
//      r3.za = Qa[4];
//      r3.zv = Qv[3] + r3.za*(t - ts[3]);
//      r3.z0 = Qs[3] + Qv[3]*(t - ts[3]) + 0.5*r3.za*SQ(t - ts[3]);
//      //r3.z0 = th1[3] + th1v[3]*(t - ts[3]) + 0.5*r3.za*pow((t - ts[3]), 2);
//    }
//  if (2 == cur)
//    return;
//
//  if (ts[4] < t && t < ts[5])
//    {
//      r3.za = Qa[5];
//      r3.zv = Qv[4] + r3.za*(t - ts[4]);
//      r3.z0  = Qs[4] + Qv[4]*(t - ts[4]) + 0.5*r3.za*SQ(t - ts[4]);
//      //	r3.z0  = th1[4] + th1v[4]*(t - ts[4]) + 0.5*r3.za*pow((t - ts[4]), 2);
//    }
//  if (3 == cur)
//    return;
//
//  if (ts[5] < t && t < ts[6])
//    {
//      r3.za = Qa[6];
//      r3.zv = Qv[5] + r3.za*(t - ts[5]);
//      r3.z0 = Qs[5] + Qv[5]*(t - ts[5]) + 0.5*r3.za*SQ(t - ts[5]);
//      //r3.z0 = th1[5] + th1v[5]*(t - ts[5]) + 0.5*r3.za*pow((t - ts[5]), 2);
//    }
//  if (ts[6] < t && t < ts[7])
//    {
//      r3.za = Qa[7];
//      r3.zv = Qv[6] + r3.za*(t - ts[6]);
//      r3.z0 = Qs[6] + Qv[6]*(t - ts[6]) + 0.5*r3.za*SQ(t - ts[6]);
//      //	r3.z0 = th1[6] + th1v[6]*(t - ts[6]) + 0.5*r3.za*pow((t - ts[6]), 2);
//    }
//  if (ts[7] < t && t < ts[8])
//    {
//      r3.za = Qa[8];
//      r3.zv = Qv[7] + r3.za*(t - ts[7]);
//      r3.z0 = Qs[7] + Qv[7]*(t - ts[7]) + 0.5*r3.za*SQ(t - ts[7]);
//      //	r3.z0  = th1[7] + th1v[7]*(t - ts[7]) + 0.5*r3.za*pow((t - ts[7]), 2);
//    }
//  if (ts[8] <= t)
//    {// Through path
//    }
//  if (4 == cur)
//    return;
//#undef SQ
//}

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
		  //		  Mat[current_step][na2][na3] = na1;
		}
	      if (enemin2 > Energ[na1][na2][na3])
		{
		  enemin2 = Energ[na1][na2][na3];
		  //		  sat[6] = na3; //Magick max_cur(= n_divide-2-2)
		  //		  sat[5] = na2; //Magick max_cur-1
		  //		  sat[4] = na1; //Magick max_cur
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

double IDP::Hikaku2(int current_step)
{
  //assert(current_step <= (n_divide-2)-2);
  //int na1, na2, na3, na4, na5, na6;
  int nx1, nx2, nx3, ny1, ny2, ny3;
  enemin2 = DBL_MAX;
  
  for(nx1 = 1; nx1 <= NN; nx1++){
    for(ny1 = 1; ny1 <= NN; ny1++){
      // --------------------------------------------------
     for(nx2 = 1; nx2 <= NN; nx2++){
	for(ny2 = 1; ny2 <= NN; ny2++){
	  // --------------------------------------------------
	  enemin1 = DBL_MAX;
	  for(nx3 = 1; nx3 <= NN; nx3++){
	    for(ny3 = 1; ny3 <= NN; ny3++){
	      // --------------------------------------------------
	      if (enemin1 > Energ2[nx1][ny1][nx2][ny2][nx3][ny3])
		{
		  enemin1 = Energ2[nx1][ny1][nx2][ny2][nx3][ny3];
		  //Mat[current_step][na2][na3] = na1;
		  Mxt[current_step][nx2][ny2][nx3][ny3] = nx1;
		  Myt[current_step][nx2][ny2][nx3][ny3] = ny1;
		}
	      if (enemin2 > Energ2[nx1][ny1][nx2][ny2][nx3][ny3])
		{
		  enemin2 = Energ2[nx1][ny1][nx2][ny2][nx3][ny3];
		  // sxt[6] = nx3; sxt[5] = nx2; sxt[4] = nx1;
		  // sxt[6] = nx3; sxt[5] = nx2; sxt[4] = nx1;
		  sxt[6] = nx3; sxt[5] = nx2; sxt[4] = nx1;
		  syt[6] = ny3; syt[5] = ny2; syt[4] = ny1;
		  // szt[6] = na3; szt[5] = na2; szt[4] = na1;
		  //sat[n_divide - 2] = na3; //Magick max_cur(= n_divide-2-2)
		  //sat[n_divide - 1] = na2; //Magick max_cur-1
		  //sat[n_divide - 0] = na1; //Magick max_cur
		  //sat[6] = na3; //Magick max_cur(= n_divide-2-2)
		  //sat[5] = na2; //Magick max_cur-1
		  //sat[4] = na1; //Magick max_cur
		}
	    } // step ny3
	  }   // step ny2
	}     // step ny1
      }	      // step nx3
    }	      // step nx2
  }	      // step nx1
  //assert(enemin1 >= 0 && enemin2 >= 0);
  return enemin2;
}

// double IDP::Hikaku3(int current_step)
// {
//   //assert(current_step <= (n_divide-2)-2);
//   //int na1, na2, na3, na4, na5, na6;
//   int nx1, nx2, nx3, ny1, ny2, ny3, nz1, nz2, nz3;
//   enemin2 = DBL_MAX;
  
//   for(nx1 = 1; nx1 <= NN; nx1++){
//     for(ny1 = 1; ny1 <= NN; ny1++){
//       for(nz1 = 1; nz1 <= NN; nz1++){
// 	// --------------------------------------------------
// 	for(nx2 = 1; nx2 <= NN; nx2++){
// 	  for(ny2 = 1; ny2 <= NN; ny2++){
// 	    for(nz2 = 1; nz2 <= NN; nz2++){
// 	      // --------------------------------------------------	  
// 	      enemin1 = DBL_MAX;		
// 	      for(nx3 = 1; nx3 <= NN; nx3++){
// 		for(ny3 = 1; ny3 <= NN; ny3++){
// 		  for(nz3 = 1; nz3 <= NN; nz3++){
// 		    // --------------------------------------------------	      
// 		    if (enemin1 > Energ3[nx1][ny1][nz1][nx2][ny2][nz2][nx3][ny3][nz3])
// 		      {
// 			enemin1 = Energ3[nx1][ny1][nz1][nx2][ny2][nz2][nx3][ny3][nz3];
// 			//Mat[current_step][na2][na3] = na1;
// 			Mxt[ current_step ][ nx2 ][ nx3 ] = nx1;
// 			Myt[ current_step ][ ny2 ][ ny3 ] = ny1;
// 			Mzt[ current_step ][ nz2 ][ nz3 ] = nz1;
// 		      }
// 		    if (enemin2 > Energ3[nx1][ny1][nz1][nx2][ny2][nz2][nx3][ny3][nz3])
// 		      {
// 			enemin2 = Energ[nx1][ny1][nz1][nx2][ny2][nz2][nx3][ny3][nz3];
// 			sxt[6] = nx3; sxt[5] = nx2; sxt[4] = nx1;
// 			syt[6] = ny3; syt[5] = ny2; syt[4] = ny1;
// 			szt[6] = nz3; szt[5] = nz2; szt[4] = nz1;
// 		      }
// 		  } // nz3
// 		}   // ny3
// 	      }	    // nx3
// 	      // --------------------------------------------------
// 	    } // nz2
// 	  }   // ny2
// 	}     // nx2
// 	// --------------------------------------------------
//       }	// nz1
//     }	// ny1
//   }	// nx1

//   //assert(enemin1 >= 0 && enemin2 >= 0);
//   return enemin2;
// }

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
//   DISPLAY(os, Qs_i);	// not initialize
//   DISPLAY(os, Qs_f);
//   DISPLAY(os, Qv_i);	// not initialize
//   DISPLAY(os, Qv_f);
//   DISPLAY(os, Qa_i);	// not initialize
//   DISPLAY(os, Qa_f);    
  // Program terminated with SIGFPE,
  // Arithmetic exception.
  //  ARYDISPLAY(os, Sth );
  //ARYDISPLAY(os, th );
  ARYDISPLAY(os, ntmax );
  DISPLAY(os, NN);
  DISPLAY(os, N3);
//   ARYDISPLAY(os, Qs);
//   ARYDISPLAY(os, Qv);
//   ARYDISPLAY(os, Qa);    
//   DISPLAY(os, kiza);
  //DISPLAY(os, t);
}
#undef DISPLAY
#undef ARYDISPLAY

const Model& IDP::getState() const
{
  return state_w;
}

//double IDPUP::En()
//{
//  e3.tau = r3.m*(r3.za + g)*r3.radius;// pulley's radius and torque.
//  e3.ev  = b1*(r3.zv/r3.radius) + b2*(r3.za/r3.radius) + b3*e3.tau;
//  // Nov18 19:00
//  // follow paragraph is comment out 
//  // Nov19
//  //if (r3.zv > 0.0) {
//  //  e3.ev = e3.ev + b3*Mf;
//  //} else if (r3.zv < 0.0) {
//  //  e3.ev = e3.ev - b3*Mf;
//  //}
//  e3.ia  = (e3.ev - Kv*r3.zv/r3.radius)/Ra;
//  
//  if (e3.ev*e3.ia > 0.0){
//    e3.energy = e3.energy + e3.ev*e3.ia*dt;
//  }
//  return e3.energy;
//}

void
IDP::set_position(){}

//IDPCurve::IDPCurve(){}
//IDPCurve::~IDPCurve(){}
//IDPCurve::IDPCurve(string _name,
//		   double _t_i, double _t_f,
//		   double _dt, double _kiza,
//		   int _n_repeat, int _n_divide, int _NN)
//  : IDP( _name,	 _t_i,  _t_f, _dt,  _kiza, _n_repeat, _n_divide, _NN)
//{
//  debug();
//  sprintf(fname_r1_sinc_st,
//	  "%s/%s_r1_sinc_st.dat"
//	  , dname_dir_etc, name.c_str() );
//  sprintf(fname_r2_sinc_st,
//	  "%s/%s_r2_sinc_st.dat"
//	  , dname_dir_etc, name.c_str() );
//  sprintf(fname_r1_sinc_ev,
//	  "%s/%s_r1_sinc_ev.dat"
//	  , dname_dir_etc, name.c_str() );
//  sprintf(fname_r2_sinc_ev,
//	  "%s/%s_r2_sinc_ev.dat"
//	  , dname_dir_etc, name.c_str() );
//
//  sprintf(ftemp_r1_saiteki_st, "%s/%s_r1_saiteki_st_%%05d.dat"
//	  , dname_dir_etc, name.c_str() );
//  sprintf(ftemp_r2_saiteki_st, "%s/%s_r2_saiteki_st_%%05d.dat"
//	  , dname_dir_etc, name.c_str() );
//  sprintf(ftemp_r1_saiteki_ev, "%s/%s_r1_saiteki_ev_%%05d.dat"
//	  , dname_dir_etc, name.c_str() );
//  sprintf(ftemp_r2_saiteki_ev, "%s/%s_r2_saiteki_ev_%%05d.dat"
//	  , dname_dir_etc, name.c_str() );
//
//  sprintf(fname_r1_enemin2, "%s/%s_r1_enemin2.dat"
//	  , dname_dir_etc, name.c_str() );
//  sprintf(fname_r2_enemin2, "%s/%s_r2_enemin2.dat"
//	  , dname_dir_etc, name.c_str() );
//}

//IDPCurve::Init(Model _state_i,
void
IDP::Init(Model _state_i,
	  SinCurv _sc, Bspline2d _bs,
	  Point pos_i, Point pos_f)
{
  printf("%s-%d\n", __FUNCTION__, __LINE__);
  
  state_i = _state_i;
  state_w = _state_i;
  sc = _sc;
  bs = _bs;
  
  Xs_i = pos_i.x; Xs_f = pos_f.x;
  Ys_i = pos_i.y; Ys_f = pos_f.y;
  Zs_i = pos_i.z; Zs_f = pos_f.z;
  

  double para = 0.0;
  //  for(int i = 0; i <= n_divide - 2; ++i)
  for(int i = 0; i <= n_divide ; ++i)  
    {
      para = sc.sin_th1( ts[i] );
      Sxx[i] = bs.X( para );
      Syy[i] = bs.Y( para );
      Szz[i] = bs.Z( para );
    }
  //FILE *fp = fopen("dir-etc/kiseki.dat","w");
  //for(int i = 0; i <= n_divide; ++i)
  //  {
  //    fprintf(fp, "%lf %lf %lf %lf %lf\n"
  //	      , Sxx[i], Syy[i], Szz[i], t1*i, ts[i]);
  //  }
  //fclose(fp);
  
  //assert(false);
}

void IDP::sinc()
{
  char fname_energy[64]; sprintf(fname_energy, ftemp_energy, 0);
  ofstream ofs_energy(fname_energy);
  ofstream ofs_r1_st(fname_r1_sinc_st), ofs_r1_ev(fname_r1_sinc_ev);
  ofstream ofs_r2_st(fname_r2_sinc_st), ofs_r2_ev(fname_r2_sinc_ev);
  
  double t = 0;
  double para[3]={0.0, 0.0, 0.0};
  double vxs[3], vys[3], vzs[3];
  double vxv[2], vyv[2], vzv[2];
  double vxa[1], vya[1], vza[1];
  
  double step = 1.0/static_cast<int>(Ts/dt);
  for( int i = 0; i <= static_cast<int>(Ts/dt) - 2; ++i )
    {
      // real manipulator 1
      graph_standard(ofs_r1_st, t,
		     r1.th, r1.thv, r1.tha,
		     e1.ev, e1.ia,  e1.tau,
		     e1.energy, e1.ev*e1.ia);
      graph_volt_constituent(ofs_r1_ev, t, e1.ev,
			     b1*r1.thv, b2*r1.tha, b3*e1.tau);
      
      // real manipulator 2
      graph_standard(ofs_r2_st, t,
		     r2.th, r2.thv, r2.tha,
		     e2.ev, e2.ia, e2.tau,
		     e2.energy, e2.ev*e2.ia);
      graph_volt_constituent(ofs_r2_ev, t, e2.ev,
			     b1*r2.thv, b2*r2.tha, b3*e2.tau);
      
      // real manipulator 3

      t	= t + dt;
      
      para[0] = para[0] + step;
      para[1] = para[0] + step;
      para[2] = para[1] + step;      
      // position--------------------------------------------------
      vxs[2] = bs.X( para[2] );
      vxs[1] = bs.X( para[1] );
      vxs[0] = bs.X( para[0] );
      
      vys[2] = bs.Y( para[2] );      
      vys[1] = bs.Y( para[1] );
      vys[0] = bs.Y( para[0] );
      
      vzs[2] = bs.Z( para[2] );
      vzs[1] = bs.Z( para[1] );
      vzs[0] = bs.Z( para[0] );
      // velocity--------------------------------------------------      
      vxv[1] = (vxs[2] - vxs[1])/dt;
      vyv[1] = (vys[2] - vys[1])/dt;
      vzv[1] = (vzs[2] - vzs[1])/dt;
      
      vxv[0] = (vxs[1] - vxs[0])/dt;
      vyv[0] = (vys[1] - vys[0])/dt;
      vzv[0] = (vzs[1] - vzs[0])/dt;
      // acceleration--------------------------------------------------
      vxa[0] = (vxv[1] - vxv[0])/dt;
      vya[0] = (vyv[1] - vyv[0])/dt;
      vza[0] = (vzv[1] - vzv[0])/dt;
      
      v1.x = vxs[0]; v1.xv = vxv[0]; v1.xa = vxa[0];
      v1.y = vys[0]; v1.yv = vyv[0]; v1.ya = vya[0];
      v1.z = vzs[0]; v1.zv = vzv[0]; v1.za = vza[0];
      
      enemin2 = En();
      
      ofs_energy << t <<" "<< enemin2 << endl;
      
      if (i % 10 == 0)
	{
	  char fname[4][64];
	  sprintf(fname[1], ftemp_sinc_e[1], i);
	  sprintf(fname[2], ftemp_sinc_e[2], i);	
	  sprintf(fname[3], ftemp_sinc_e[3], name.c_str(), i);
	  ofstream ofs1(fname[1]);
	  ofstream ofs2(fname[2]);
	  ofstream ofs3(fname[3]);
	  //	  zu( ofs1 ); //All manipulator draw red color.
	  zu( ofs1, ofs2 );
	  // zu( ofs1, ofs2, ofs3 );// !not implement
	}
      // output volt ampere field
      // plot("");
    }

  for(int i = 1; i <= 2; ++i)
    {
      // real manipulator 1
      graph_standard(ofs_r1_st, t,
		     r1.th, r1.thv, r1.tha,
		     e1.ev, e1.ia,  e1.tau,
		     e1.energy, e1.ev*e1.ia);
      graph_volt_constituent(ofs_r1_ev, t, e1.ev,
			     b1*r1.thv, b2*r1.tha, b3*e1.tau);
      
      // real manipulator 2
      graph_standard(ofs_r2_st, t,
		     r2.th, r2.thv, r2.tha,
		     e2.ev, e2.ia, e2.tau,
		     e2.energy, e2.ev*e2.ia);
      graph_volt_constituent(ofs_r2_ev, t, e2.ev,
			     b1*r2.thv, b2*r2.tha, b3*e2.tau);

      t = t + dt;
      v1.x = vxs[i]; v1.xv = vxv[1]; v1.xa = vxa[0];
      v1.y = vys[i]; v1.yv = vyv[1]; v1.yv = vya[0];
      v1.z = vzs[i]; v1.zv = vzv[1]; v1.za = vza[0];
      
      enemin2 = En();
      ofs_energy << t <<" "<< enemin2 << endl;
    }
  //assert(false);
}

//void IDPCurve::settei(int cur)
void IDP::settei(int cur)
{//  assert(cur <= (n_divide-2)-2);
  Xv[1] = 2.0*(Xs[1] - Xs_i)/t1;
  Xa[1] = Xv[1]/t1;
  Xv[2] = 2.0*(Xs[2] - Xs[1])/t1 - Xv[1];
  Xa[2] = (Xv[2] - Xv[1])/t1;
  Xv[3] = 2.0*(Xs[3] - Xs[2])/t1 - Xv[2];
  Xa[3] = (Xv[3] - Xv[2])/t1;

  Yv[1] = 2.0*(Ys[1] - Ys_i)/t1;
  Ya[1] = Yv[1]/t1;
  Yv[2] = 2.0*(Ys[2] - Ys[1])/t1 - Yv[1];
  Ya[2] = (Yv[2] - Yv[1])/t1;
  Yv[3] = 2.0*(Ys[3] - Ys[2])/t1 - Yv[2];
  Ya[3] = (Yv[3] - Yv[2])/t1;

  //   Zv[1] = 2.0*(Zs[1] - Zs_i)/t1;
  //   Za[1] = Zv[1]/t1;
  //   Zv[2] = 2.0*(Zs[2] - Zs[1])/t1 - Zv[1];
  //   Za[2] = (Zv[2] - Zv[1])/t1;
  //   Zv[3] = 2.0*(Zs[3] - Zs[2])/t1 - Zv[2];
  //   Za[3] = (Zv[3] - Zv[2])/t1;
  if (cur == 1) return;
    
  Xv[4] = 2.0*(Xs[4] - Xs[3])/t1 - Xv[3];
  Xa[4] = (Xv[4] - Xv[3])/t1;

  Yv[4] = 2.0*(Ys[4] - Ys[3])/t1 - Yv[3];
  Ya[4] = (Yv[4] - Yv[3])/t1;

  //   Zv[4] = 2.0*(Zs[4] - Zs[3])/t1 - Zv[3];
  //   Za[4] = (Zv[4] - Zv[3])/t1;
  if (cur == 2) return;
    
  Xv[5] = 2.0*(Xs[5] - Xs[4])/t1 - Xv[4];
  Xa[5] = (Xv[5] - Xv[4])/t1;

  Yv[5] = 2.0*(Ys[5] - Ys[4])/t1 - Yv[4];
  Ya[5] = (Yv[5] - Yv[4])/t1;

  //   Zv[5] = 2.0*(Zs[5] - Zs[4])/t1 - Zv[4];
  //   Za[5] = (Zv[5] - Zv[4])/t1;
  if (cur == 3) return;
    
  Xv[6] = 2.0*(Xs[6] - Xs[5])/t1 - Xv[5];
  Xa[6] = (Xv[6] - Xv[5])/t1;

  Yv[6] = 2.0*(Ys[6] - Ys[5])/t1 - Yv[5];
  Ya[6] = (Yv[6] - Yv[5])/t1;

  //   Zv[6] = 2.0*(Zs[6] - Zs[5])/t1 - Zv[5];
  //   Za[6] = (Zv[6] - Zv[5])/t1;
  // will become if (step == 4) return;
  
  Xs[8] = Xs_f;
  Xv[8] = 0.0; // vary with condition
  Xv[7] = (Xs[8] - Xs[6])/t1 - (Xv[6]+Xv[8])/2.0;
  Xa[7] = (Xv[7] - Xv[6])/t1;
  Xs[7] = Xs[6] + Xv[6] * t1 + Xa[7]*t1*t1/2.0;
  Xa[8] = (Xv[8] - Xv[7])/t1;

  Ys[8] = Ys_f;
  Yv[8] = 0.0; // vary with condition
  Yv[7] = (Ys[8] - Ys[6])/t1 - (Yv[6]+Yv[8])/2.0;
  Ya[7] = (Yv[7] - Yv[6])/t1;
  Ys[7] = Ys[6] + Yv[6] * t1 + Ya[7]*t1*t1/2.0;
  Ya[8] = (Yv[8] - Yv[7])/t1;

  //   Zs[8] = Zs_f;
  //   Zv[8] = 0.0; // vary with condition
  //   Zv[7] = (Zs[8] - Zs[6])/t1 - (Zv[6]+Zv[8])/2.0;
  //   Za[7] = (Zv[7] - Zv[6])/t1;
  //   Zs[7] = Zs[6] + Zv[6] * t1 + Za[7]*t1*t1/2.0;
  //   Za[8] = (Zv[8] - Zv[7])/t1;
  if (cur == 4) return;
}

// if (cur == 5?)
//  return; ((n_divide-2)-3) + 1 < mat_max
// if (cur <= mat_max)
//     settei1()
// else
//     settei2(get the last two step acceration and velocity)
// 
// mat_max = ((n_divide - 2) - 3) + 1 (n_divide=8; 4);

//IDPCurve::sim(int cur, double t)

void IDP::sim(int cur, double t)
{//  assert(cur <= (n_divide-2)-2);
#define SQ(x) ((x)*(x))
  if (t < ts[1])
    {
      v1.xa = Xa[1];
      v1.xv = v1.xa*t;
      v1.x = 0.5*v1.xa*t*t + Xs_i;

      v1.ya = Ya[1];
      v1.yv = v1.ya*t;
      v1.y = 0.5*v1.ya*t*t + Ys_i;
      
      //       v1.za = Za[1];
      //       v1.zv = v1.za*t;
      //       v1.z = 0.5*v1.za*t*t + Zs_i;
    }
  
  if (ts[1] < t && t < ts[2])
    {
      v1.xa = Xa[2];
      v1.xv = Xv[1] + v1.xa*(t - ts[1]);
      v1.x  = Xs[1] + Xv[1]*(t - ts[1]) + 0.5*v1.xa*SQ(t - ts[1]);

      v1.ya = Ya[2];
      v1.yv = Yv[1] + v1.ya*(t - ts[1]);
      v1.y  = Ys[1] + Yv[1]*(t - ts[1]) + 0.5*v1.ya*SQ(t - ts[1]);

      //       v1.za = Za[2];
      //       v1.zv = Zv[1] + v1.za*(t - ts[1]);
      //       v1.z = Zs[1] + Zv[1]*(t - ts[1]) + 0.5*v1.za*SQ(t - ts[1]);
    }
  
  if (ts[2] < t && t < ts[3])
    {
      v1.xa = Xa[3];
      v1.xv = Xv[2] + v1.xa*(t - ts[2]);
      v1.x  = Xs[2] + Xv[2]*(t - ts[2]) + 0.5*v1.xa*SQ(t - ts[2]);

      v1.ya = Ya[3];
      v1.yv = Yv[2] + v1.ya*(t - ts[2]);
      v1.y  = Ys[2] + Yv[2]*(t - ts[2]) + 0.5*v1.ya*SQ(t - ts[2]);

      //       v1.za = Za[3];
      //       v1.zv = Zv[2] + v1.za*(t - ts[2]);
      //       v1.z = Zs[2] + Zv[2]*(t - ts[2]) + 0.5*v1.za*SQ(t - ts[2]);	
    }
  if (1 == cur)
    return;

  if (ts[3] < t && t < ts[4])
    {
      v1.xa = Xa[4];
      v1.xv = Xv[3] + v1.xa*(t - ts[3]);
      v1.x  = Xs[3] + Xv[3]*(t - ts[3]) + 0.5*v1.xa*SQ(t - ts[3]);

      v1.ya = Ya[4];
      v1.yv = Yv[3] + v1.ya*(t - ts[3]);
      v1.y  = Ys[3] + Yv[3]*(t - ts[3]) + 0.5*v1.ya*SQ(t - ts[3]);

      //       v1.za = Za[4];
      //       v1.zv = Zv[3] + v1.za*(t - ts[3]);
      //       v1.z = Zs[3] + Zv[3]*(t - ts[3]) + 0.5*v1.za*SQ(t - ts[3]);
    }
  if (2 == cur)
    return;

  if (ts[4] < t && t < ts[5])
    {
      v1.xa = Xa[5];
      v1.xv = Xv[4] + v1.xa*(t - ts[4]);
      v1.x  = Xs[4] + Xv[4]*(t - ts[4]) + 0.5*v1.xa*SQ(t - ts[4]);

      v1.ya = Ya[5];
      v1.yv = Yv[4] + v1.ya*(t - ts[4]);
      v1.y  = Ys[4] + Yv[4]*(t - ts[4]) + 0.5*v1.ya*SQ(t - ts[4]);

      //       v1.za = Za[5];
      //       v1.zv = Zv[4] + v1.za*(t - ts[4]);
      //       v1.z  = Zs[4] + Zv[4]*(t - ts[4]) + 0.5*v1.za*SQ(t - ts[4]);
    }
  if (3 == cur)
    return;

  if (ts[5] < t && t < ts[6])
    {
      v1.xa = Xa[6];
      v1.xv = Xv[5] + v1.xa*(t - ts[5]);
      v1.x  = Xs[5] + Xv[5]*(t - ts[5]) + 0.5*v1.xa*SQ(t - ts[5]);

      v1.ya = Ya[6];
      v1.yv = Yv[5] + v1.ya*(t - ts[5]);
      v1.y  = Ys[5] + Yv[5]*(t - ts[5]) + 0.5*v1.ya*SQ(t - ts[5]);

      //       v1.za = Za[6];
      //       v1.zv = Zv[5] + v1.za*(t - ts[5]);
      //       v1.z = Zs[5] + Zv[5]*(t - ts[5]) + 0.5*v1.za*SQ(t - ts[5]);
    }
  
  if (ts[6] < t && t < ts[7])
    {
      v1.xa = Xa[7];
      v1.xv = Xv[6] + v1.xa*(t - ts[6]);
      v1.x  = Xs[6] + Xv[6]*(t - ts[6]) + 0.5*v1.xa*SQ(t - ts[6]);

      v1.ya = Ya[7];
      v1.yv = Yv[6] + v1.ya*(t - ts[6]);
      v1.y  = Ys[6] + Yv[6]*(t - ts[6]) + 0.5*v1.ya*SQ(t - ts[6]);

      // v1.za = Za[7];
      // v1.zv = Zv[6] + v1.za*(t - ts[6]);
      // v1.z = Zs[6] + Zv[6]*(t - ts[6]) + 0.5*v1.za*SQ(t - ts[6]);
    }
  
  if (ts[7] < t && t < ts[8])
    {
      v1.xa = Xa[8];
      v1.xv = Xv[7] + v1.xa*(t - ts[7]);
      v1.x  = Xs[7] + Xv[7]*(t - ts[7]) + 0.5*v1.xa*SQ(t - ts[7]);

      v1.ya = Ya[8];
      v1.yv = Yv[7] + v1.ya*(t - ts[7]);
      v1.y  = Ys[7] + Yv[7]*(t - ts[7]) + 0.5*v1.ya*SQ(t - ts[7]);

      // v1.za = Za[8];
      // v1.zv = Zv[7] + v1.za*(t - ts[7]);
      // v1.z = Zs[7] + Zv[7]*(t - ts[7]) + 0.5*v1.za*SQ(t - ts[7]);
    }
  
  if (ts[8] <= t)
    {// Through path
    }
  if (4 == cur)
    return;
#undef SQ
}

//double IDPCurve::det(double a11, double a12, double a21, double a22)
double IDP::det(double a11, double a12, double a21, double a22)
{
  return a11*a22 - a12*a21;
}

//IDPCurve::En()
double
IDP::En()
{
  double detJ, J11, J12, J21, J22;//, coef;
  double cth1, cth12, cth2;
  double sth1, sth12, sth2;  
  double keirosokudo, PL3, XXa, YYa;

  cth1 = cos(r1.th); cth12 = cos(r1.th + r2.th); cth2 = cos(r2.th);
  sth1 = sin(r1.th); sth12 = sin(r1.th + r2.th); sth2 = sin(r2.th);

  keirosokudo = sqrt(v1.xv*v1.xv + v1.yv*v1.yv);
  J11 = -r1.l*sth1 - r2.l*sth12;
  J12 = -r2.l*sth12;
  J21 =  r1.l*cth1 + r2.l*cth12;
  J22 =  r2.l*cth12;
  
  detJ = det(J11, J12, J21, J22);
  r1.thv = det(v1.xv, J12 ,
	       v1.yv, J22 )/detJ;// => arithmetic exception
  r2.thv = det(J11, v1.xv,
	       J21, v1.yv)/detJ;

  PL3 = sqrt(v1.x*v1.x + v1.y*v1.y); // v1.x ? (v1.x - v1.x0)
  r1.th = atan2(v1.y, v1.x) - acos(0.5*PL3/r1.l);
  r2.th = 2.0*acos(0.5*PL3/r1.l);

  XXa = v1.xa + r1.thv*r1.thv*r1.l*cth1
    + (r1.thv + r2.thv)*(r1.thv + r2.thv)*r2.l*cth12;
  YYa = v1.ya + r1.thv*r1.thv*r1.l*sth1
    + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*sth12;

  r1.tha = det(XXa, J12,
	       YYa, J22)/detJ;
  r2.tha = det(J11, XXa,
	       J21, YYa)/detJ;
  
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
    
#define WATT( e )  ( (e).ev*(e).ia )
  if ( WATT(e1) > 0.0)
    {
      e1.energy = e1.energy + WATT(e1)*dt;
    }
  if ( WATT(e2) > 0.0)
    { //      Ene_2 = Ene_2 + ev2*ia2*dt;
      e2.energy = e2.energy + WATT(e2)*dt;
    }
  // if ( WATT(e3) > 0.0){ //      Ene_3 = Ene_3 + ev3*ia3*dt;
  //   e3.energy = e3.energy + WATT(e3)*dt;
  // }
#undef WATT
  //Ene123 = e1.energy + e2.energy + e3.energy;
  //  return Ene123; // => Nov 23
  //Ene123 = e1.energy + e2.energy;
  return e1.energy + e2.energy;
}

//void IDPCurve::zu(ostream& ofs, vector<Point3D>trjc )
void IDP::zu(ostream& ofs, vector<Point> trjc)
{
  GraphicalObject* obj[] =
    {
      new Circle3D(v1.x0, v1.y0, v1.z0, 0.02),
      new Circle3D(v1.x , v1.y , v1.z0 , 0.01),
      new Line3D(v1.x0, v1.y0, v1.z0, v1.x , v1.y , v1.z0 ),
    };

  for(uint i = 0; i < sizeof(obj)/sizeof(*obj); ++i)
    {
      obj[i]->draw(ofs);
    }
  
  for(uint i = 0; i < trjc.size(); ++i)
    {
      ofs << trjc[i].x <<"\t"
	  << trjc[i].y <<"\t"
	  << trjc[i].z << endl;
    }
}

//void IDPCurve::zu(ostream& ofs)
void IDP::zu(ostream& ofs)
{
  static double rr = (pnt_f.x - pnt_i.x)/10.0;
  r1.x = r1.x0 + r1.l*cos(r1.th);
  r1.y = r1.y0 + r1.l*sin(r1.th);
  
  r2.x0 = r1.x;
  r2.y0 = r1.y;
  
  r2.x = r2.x0 + r2.l*cos(r1.th + r2.th);
  r2.y = r2.y0 + r2.l*sin(r1.th + r2.th);
  
  r3.x0 = r2.x;
  r3.y0 = r2.y;
  r3.z0 = r3.z + r3.l;
  
  GraphicalObject* obj[] =
    {
      new Circle3D(r1.x0, r1.y0, r1.z0, rr), // 3
      new Circle3D(r1.x , r1.y , r1.z0, rr), // 4
      new Circle3D(r2.x0, r2.y0, r2.z0, rr), // 4
      new Circle3D(r3.x0, r3.y0, r3.z0, rr), // 5
      new Circle3D(r3.x0, r3.y0, r3.z , rr), // 6
      new Line3D(r1.x0, r1.y0, r1.z0, r1.x, r1.y, r1.z0),   // 7
      new Line3D(r2.x0, r2.y0, r2.z0, r2.x, r2.y, r2.z0),   // 8
      new Line3D(r3.x0, r3.y0, r3.z0, r3.x0, r3.y0, r3.z ), // 9
    };

  for(uint i = 0; i < sizeof(obj)/sizeof(*obj); ++i)
    {
      obj[i]->draw(ofs);
    }
}

//void IDPCurve::zu(ostream& ofs1, ostream& ofs2)
void IDP::zu(ostream& ofs1, ostream& ofs2)
{
  static double rr = (pnt_f.x - pnt_i.x)/10.0;

  r1.x = r1.x0 + r1.l*cos(r1.th);
  r1.y = r1.y0 + r1.l*sin(r1.th);
  
  r2.x0 = r1.x; r2.x = r2.x0 + r2.l*cos(r1.th + r2.th);
  r2.y0 = r1.y; r2.y = r2.y0 + r2.l*sin(r1.th + r2.th);
  
  r3.x0 = r2.x;
  r3.y0 = r2.y;
  r3.z0 = r3.z - r3.l;
  
  GraphicalObject* obj[] =
    {
      new Circle3D(v1.x0, v1.y0, v1.z0, rr*0.8), // 0
      new Circle3D(v1.x , v1.y , v1.z , rr*0.5), // 1
      new Line3D(v1.x0, v1.y0, v1.z0, v1.x , v1.y , v1.z), // 2
      // -------------------------------------------------- //
      new Circle3D(r1.x0, r1.y0, r1.z0, rr), // 3
      new Circle3D(r2.x0, r2.y0, r2.z0, rr), // 4
      new Circle3D(r3.x0, r3.y0, r3.z0, rr), // 5
      new Circle3D(r3.x0, r3.y0, r3.z0, rr), // 6
      new Line3D(r1.x0, r1.y0, r1.z0, r1.x, r1.y, r1.z0), // 7
      new Line3D(r2.x0, r2.y0, r2.z0, r2.x, r2.y, r2.z0), // 8
      new Line3D(r3.x0, r3.y0, r3.z0, r3.x0, r3.y0, r3.z ), // 9
    };
  
  for(uint i = 0; i < sizeof(obj)/sizeof(*obj); ++i)
    {
      if (i < 3)
	{
	  obj[i]->draw(ofs1);
	}
      else
	{
	  obj[i]->draw(ofs2);
	}
    }
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

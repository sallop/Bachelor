#ifndef __INCLUDE_IDP_H__
#define __INCLUDE_IDP_H__
#include <string>
#include <cstring>
#include <vector>
#include "GraphicalPrimitive.hpp"
#include "Manipulator.hpp"
#include "Energy.hpp"
#include "variables.hpp"
#include "Bi3_2.hpp"

using namespace std;

class IDP;
//class IDPCurve;
class Model;

struct Point3D
{
  double x, y, z;
  Point3D(double _x, double _y, double _z): x(_x), y(_y), z(_z){}
};

class Model
{
public:
  Model();
  Model(VRManipulator _v1,
	Revolute _r1, Revolute _r2, Prismatic _r3,
	Energy _e1, Energy _e2, Energy _e3);
  ~Model();
  VRManipulator v1;
  Revolute	r1;
  Revolute	r2;
  Prismatic	r3;
  Energy	e1;
  Energy	e2;
  Energy	e3;
};

class IDP
{
protected:
  int make_animation_flag[60];
  char dname_dir_dat[64];
  char dname_dir_etc[64];
  
  char fname_sinc_st[64];
  char fname_sinc_ev[64];

  char ftemp_saiteki_st[64];
  char ftemp_saiteki_ev[64];

  char fname_enemin2[64];
  char ftemp_energy[64];
  
  char ftemp_sinc_e[4][64];
  char ftemp_saiteki_e[4][64];

  char fname_r1_sinc_st[64];
  char fname_r2_sinc_st[64];
  char fname_r3_sinc_st[64];
  char fname_r1_sinc_ev[64];
  char fname_r2_sinc_ev[64];
  char fname_r3_sinc_ev[64];
  char ftemp_r1_saiteki_st[64];
  char ftemp_r2_saiteki_st[64];
  char ftemp_r3_saiteki_st[64];  
  char ftemp_r1_saiteki_ev[64];
  char ftemp_r2_saiteki_ev[64];
  char ftemp_r3_saiteki_ev[64];  
  char fname_r1_enemin2[64];
  char fname_r2_enemin2[64];
  char fname_r3_enemin2[64];
  
  //static int arysz;		// default array size
  enum{arysz = 20};
  Model state_i, state_w;
  
  VRManipulator &v1;		// alias v1_w
  Revolute	&r1, &r2;	// alias state_w.r1
  Prismatic&	 r3;		// alias state_w.r3
  Energy	&e1, &e2, &e3;	// alias state_w.e1

#define DIVIDE	8
#define NNN	7
#define SATSZ	6
#define CURRENT 4
  //  int at[arysz], sat[arysz];
  
  // IDP::IDP	 = > class variable
  //  int Mat[1+ 4][1+ 7][1+ 7];	// [n_divide-2-2 = max(cur)][NN][NN]
  //  double Energ[1+ 7][1+ 7][1+ 7];	// [NN][NN][NN]

  int xt[1+SATSZ], sxt[1+SATSZ];
  int yt[1+SATSZ], syt[1+SATSZ];
  int zt[1+SATSZ], szt[1+SATSZ];
  int Mxt[1+ CURRENT][1+ NNN][1+ NNN][1+ NNN][1+ NNN];
  int Myt[1+ CURRENT][1+ NNN][1+ NNN][1+ NNN][1+ NNN];
  // int Mxt[1+ CURRENT][1+ NNN][1+ NNN][1+ NNN][1+ NNN][1+ NNN][1+ NNN];
  // int Myt[1+ CURRENT][1+ NNN][1+ NNN][1+ NNN][1+ NNN][1+ NNN][1+ NNN];
  // int Mzt[1+ CURRENT][1+ NNN][1+ NNN][1+ NNN][1+ NNN][1+ NNN][1+ NNN];
  
  double Xs[1+DIVIDE], Xv[1+DIVIDE], Xa[1+DIVIDE], Xs_i, Xs_f;
  double Ys[1+DIVIDE], Yv[1+DIVIDE], Ya[1+DIVIDE], Ys_i, Ys_f;
  double Zs[1+DIVIDE], Zv[1+DIVIDE], Za[1+DIVIDE], Zs_i, Zs_f;
  
  //  double	Sth[arysz];  
  double Sxx[1+DIVIDE];
  double Syy[1+DIVIDE];
  double Szz[1+DIVIDE];
  // double Qs_f, Qs_i, Qs[arysz];
  // double Qv_f, Qv_i, Qv[arysz];
  // double Qa_f, Qa_i, Qa[arysz];
  double Energ[1+ NNN][1+ NNN][1+ NNN];
  double Energ2[1+ NNN][1+ NNN][1 + NNN][1+ NNN][1+ NNN][1 + NNN];
  //  Energy2[xt[c]][yt[c]][xt[c+1]][yt[c+1]][xt[c+2]][yt[c+2]];
  //int Energ3[1+ NNN][1+ NNN][1+ NNN][1+ NNN][1 + NNN][1+ NNN];
  //int Energ3[xt[c]][yt[c]][zt[c]][xt[c+1]][yt[c+1]][zt[c+1]]
  //          [xt[c+2]][yt[c+2]][zt[c+2]]
#undef NNN
#undef SATSZ
#undef CURRENT
public:  
  //  double	*Enemin;	// size is depend repeat
  double Enemin[120];
  int l_nnn;		// offset is by one.
  double enemin1, enemin2;	// using Hikaku1, so it will can delete.
  int n_repeat;
  int n_divide;	// 
  int NN, N3;	// IDP Test case(7, 4)
  string name;		// output file prefix
  
  SinCurv sc;
  Bspline2d bs;

  void sinc();
  void repeat();
  void tansaku();
  void	En_sim(int ntmax, int step);
  double det(double a11, double a12, double a21, double a22);
  double En();
  double optimal_result(int make_animation_flag);
  void set_position();
  void zu(ostream& os);
  void zu(ostream& os1, ostream& os2);
  void zu(ostream& os1, vector<Point> trjc);
  void settei(int cur);
  void	sim(int cur, double t);
  double Hikaku(int jj);
  double Hikaku2(int jj);
  double Hikaku3(int jj);

  const Model& getState() const;
  void makeGifFile() const;
  
  IDP();
  IDP(string		_name,
      double		_t_i, double _t_f,
      double		_dt,
      double _kiza_x, double _kiza_y, double _kiza_z,
      int		_n_repeat, int _n_divide, int _NN);
  virtual ~IDP();
  void Init(Model _state_i,
	    SinCurv _sc, Bspline2d   _bs,
	    Point pos_i, Point pos_f);
  
  double t_i, t_f, Ts, dt, t1;	// IDP
  double ts[arysz];	// offset is zero.
  int	 ntmax[arysz];	// IDP
  //  double kiza;		// IDP
  double kiza_x, kiza_y, kiza_z;

  //virtual double En(double s, double v, double a);
  void	start();
  void	dump(ostream& ofs);
};


#endif // __INCLUDE_IDP_H__

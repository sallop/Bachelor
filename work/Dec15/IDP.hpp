#ifndef __INCLUDE_IDP_H__
#define __INCLUDE_IDP_H__
#include <string>
#include <cstring>
#include "GraphicalPrimitive.hpp"
#include "Manipulator.hpp"
#include "Energy.hpp"
#include "variables.hpp"

class IDP;
class IDPUP;
class IDPDIRECT;
class IDPDOWN;
class Model;

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

  //static int arysz;		// default array size
  enum{arysz = 20};
  Model state_i, state_w;
  
  //ScaraRobot state_i;//state_f	// Init()
  //VRManipulator v1_i;//, v1_w;// v1_f;	// init, work, final
  //ScaraRobot state_w;
  //  VRManipulator v1_w;		// working state implements

  VRManipulator &v1;		// alias v1_w
  Revolute	&r1, &r2;	// alias state_w.r1
  Prismatic&	 r3;		// alias state_w.r3
  Energy	&e1, &e2, &e3;	// alias state_w.e1
  
  int		 at[arysz], sat[arysz];	// IDP::IDP	 = > class variable
  int		 Mat[1+ 4][1+ 7][1+ 7];	// [n_divide-2-2 = max(cur)][NN][NN]
  double	 Energ[1+ 7][1+ 7][1+ 7];	// [NN][NN][NN]
  //  double	*Enemin;	// size is depend repeat
  double Enemin[120];
  int		 l_nnn;		// offset is by one.
  double	 enemin1, enemin2;	// using Hikaku1, so it will can delete.
  int		 n_repeat, n_divide;	// 
  int		 NN, N3;	// IDP Test case(7, 4)
  string	 name;		// output file prefix

  virtual void sinc() = 0;
  void		repeat();
  void		tansaku();
  virtual void	En_sim(int ntmax, int step);
  virtual double En() = 0;
  double	optimal_result(int make_animation_flag);
  virtual void	set_position();
  void		zu(ostream& os);
  void		zu(ostream& os1, ostream& os2);
  void		settei(int cur);
  virtual void	sim(int cur, double t);
  double	Hikaku(int jj);

public:
  //ScaraRobot state_f;
  //VRManipulator v1_f;
  //  Model			state_f;
  const Model& getState() const;
  void makeGifFile() const;
  //  static const int arysz;// = 20;
  IDP();
  IDP(string		_name,
      double		_t_i, double _t_f,
      double		_dt, double _kiza,
      int		_n_repeat, int _n_divide, int _NN);
  virtual ~IDP();
  void Init(Model	_state_i,//ScaraRobot _state_i,VRManipulator _v1_i,
	    SinCurv	_sc,
	    double	_Qs_i, double _Qs_f,
	    double	_Qv_i, double _Qv_f,
	    double	_Qa_i, double _Qa_f);

  double	t_i, t_f, Ts, dt, t1;	// IDP
  double	ts[arysz];	// offset is zero.
  double	Sth[arysz];
  int		ntmax[arysz];	// IDP
  double	kiza;		// IDP

  double	Qs_f, Qs_i, Qs[arysz];
  double	Qv_f, Qv_i, Qv[arysz];
  double	Qa_f, Qa_i, Qa[arysz];
  SinCurv	sc;		// Init
  
  //virtual double En(double s, double v, double a);
  void	start();
  void	dump(ostream& ofs);
};

class IDPUP : public IDP
{
public:
  IDPUP();
  IDPUP(string	_name,
	double	_t_i, double _t_f,
	double	_dt, double _kiza,
	int	_n_repeat, int _n_divide, int _NN);
  ~IDPUP();
  double	En();
  void		sinc();
  void		sim(int cur, double t);
  void		set_position();
};

class IDPDIRECT : public IDP
{
  Linear_curve_t fk1;		// constant
  Linear_curve_t fk2;		// variable
  inline double det(double a11, double a12, double a13, double a14);

  char fname_r1_sinc_st[64];
  char fname_r2_sinc_st[64];
  char fname_r1_sinc_ev[64];
  char fname_r2_sinc_ev[64];
  
  char ftemp_r1_saiteki_st[64];
  char ftemp_r2_saiteki_st[64];
  char ftemp_r1_saiteki_ev[64];
  char ftemp_r2_saiteki_ev[64];

  char fname_r1_enemin2[64];
  char fname_r2_enemin2[64];
public:
  IDPDIRECT();
  IDPDIRECT(string _name,
	    double _t_i, double _t_f,
	    double _dt, double _kiza,
	    int _n_repeat, int _n_divide, int _NN);  
  ~IDPDIRECT();

  void Init(//ScaraRobot _state_i,
	    //VRManipulator _v1_i,
	    Model		_state_i,
	    SinCurv		_sc,
	    Linear_curve_t	_fk1,  Linear_curve_t _fk2,
	    double		_Qs_i, double _Qs_f,
	    double		_Qv_i, double _Qv_f,
	    double		_Qa_i, double _Qa_f);

  void sinc();
  void sim(int cur, double t);
  double En();
  void set_position();
//void Init(Model		_state_i,
//	    Model		_state_f,
//	    SinCurv		_sc,
//	    Linear_curve_t	fk1,
//	    Linear_curve_t	fk2);
};

class IDPDOWN : public IDP
{
public:
  IDPDOWN();
  IDPDOWN(string _name,
	  double _t_i, double _t_f,
	  double _dt, double _kiza,
	  int _n_repeat, int _n_divide, int _NN);
  ~IDPDOWN();
  void sinc();
  void sim(int cur, double t);
  double En();
  void set_position();
};


#endif // __INCLUDE_IDP_H__

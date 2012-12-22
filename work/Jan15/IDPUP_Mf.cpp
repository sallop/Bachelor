#include "IDP.hpp"

IDPUP_Mf::IDPUP_Mf(){}
IDPUP_Mf::~IDPUP_Mf(){}
IDPUP_Mf::IDPUP_Mf(string _name,
	     double _t_i, double _t_f,
	     double _dt, double _kiza,
	     int _n_repeat, int _n_divide, int _NN)
  :IDP( _name, _t_i,  _t_f, _dt,  _kiza, _n_repeat, _n_divide, _NN)
{
  //debug();
  printf("%s %d\n", __FUNCTION__, __LINE__);
}

void IDPUP_Mf::sinc()
{
  ofstream ofs_energy(     fname_lst[eENERGY][eSinc][eName]           );
  ofstream ofs_endeffector(fname_lst[eGRAPH_ENDEFFECTOR][eSinc][eName]);
  //ofstream ofs_st(         fname_lst[eGRAPH_ST_R0][eSinc][eName]      );
  //ofstream ofs_r1_st (     fname_lst[eGRAPH_ST_R1][eSinc][eName]);
  //ofstream ofs_r2_st (     fname_lst[eGRAPH_ST_R2][eSinc][eName]);  
  //ofstream ofs_ev(         fname_lst[eGRAPH_EV_R0][eSinc][eName]      );
  //ofstream ofs_r1_ev (     fname_lst[eGRAPH_EV_R1][eSinc][eName]);
  //ofstream ofs_r2_ev (     fname_lst[eGRAPH_EV_R2][eSinc][eName]);
  ofstream ofs_r3_st (     fname_lst[eGRAPH_ST_R3][eSinc][eName]);
  ofstream ofs_r3_ev (     fname_lst[eGRAPH_EV_R3][eSinc][eName]);
  
  double t = 0;
  for(int i = 0; i <= static_cast<int>(Ts/dt); ++i)
    {
      graph_standard(ofs_r3_st, t    ,
		     r3.z0    , r3.zv, r3.za,
		     e3.ev    , e3.ia, e3.tau,
		     e3.energy, e3.ev*e3.ia);

      graph_volt_constituent(ofs_r3_ev, t, e3.ev,
			     b1*r3.Thv(), b2*r3.Tha(), b3*e3.tau);
      
      t       = t + dt;
      r3.z0   = sc.sin_th1(t);
      r3.zv   = sc.sin_th1v(t);
      r3.za   = sc.sin_th1a(t);
      enemin2 = En();
      ofs_energy << t <<" "<< enemin2 << endl;
      if (i % 10 == 0)
	{
	  sprintf(fname_lst[eZU_E1][eSinc][eName],
		  fname_lst[eZU_E1][eSinc][eTemp], i);
	  sprintf(fname_lst[eZU_E2][eSinc][eName],
		  fname_lst[eZU_E2][eSinc][eTemp], i);
	  sprintf(fname_lst[eZU_E3][eSinc][eName],
		  fname_lst[eZU_E3][eSinc][eTemp], i);
	  
	  ofstream ofs1(fname_lst[eZU_E1][eSinc][eName]);
	  ofstream ofs2(fname_lst[eZU_E2][eSinc][eName]);
	  //	  ofstream ofs3(fname_lst[eZU_E3][eSinc][eName]);
	  //zu( ofs1 );
	  zu(ofs1, ofs2);		// real, virtual
	  //zu(ofs1, ofs2, ofs3);
	}
    }
  enemin2 = 0.0;
}

void IDPUP_Mf::sim(int cur, double t)//Magick_numbers
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

double IDPUP_Mf::En()
{
  e3.tau = r3.m*(r3.za + g)*r3.radius;// pulley's radius and torque.
  e3.ev  = b1*(r3.zv/r3.radius) + b2*(r3.za/r3.radius) + b3*e3.tau;
  // Nov18 19:00
  // follow paragraph is comment out 
  // Nov19
  if (r3.zv > 0.0)
    {
      e3.ev = e3.ev + b3*Mf;
    }
  else if (r3.zv < 0.0)
    {
      e3.ev = e3.ev - b3*Mf;
    }
  e3.ia  = (e3.ev - Kv*r3.zv/r3.radius)/Ra;
  
  if (e3.ev*e3.ia > 0.0)
    {
      e3.energy = e3.energy + e3.ev*e3.ia*dt;
    }
  return e3.energy;
}

void
IDPUP_Mf::set_position()
{
  //  v1.x0 = pnt_v0.x; v1.y0 = pnt_v0.y;
  v1.z0 = r3.z0  ;
  //  v1.x  = pnt_i.x; v1.y  = pnt_i.y;
  v1.x = v1.x0 + v1.l*cos(v1.th);
  v1.y = v1.y0 + v1.l*sin(v1.th); 
  r3.z = r3.z0 + r3.l;  
}

#include "IDP.hpp"

IDPDIRECT::IDPDIRECT(){}
IDPDIRECT::~IDPDIRECT(){}
IDPDIRECT::IDPDIRECT(string _name,
		     double _t_i, double _t_f,
		     double _dt, double _kiza,
		     int _n_repeat, int _n_divide, int _NN)
  : IDP( _name,	 _t_i,  _t_f, _dt,  _kiza, _n_repeat, _n_divide, _NN)
{
  printf("%s %d\n", __FUNCTION__, __LINE__);
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

void IDPDIRECT::sinc()
{
  //char fname_energy[64]; sprintf(fname_energy, ftemp_energy, 0);
  ofstream ofs_energy(fname_lst[eENERGY][eSinc][eName]     );
  ofstream ofs_r1_st (fname_lst[eGRAPH_ST_R1][eSinc][eName]);
  ofstream ofs_r1_ev (fname_lst[eGRAPH_EV_R1][eSinc][eName]);
  ofstream ofs_r2_st (fname_lst[eGRAPH_ST_R2][eSinc][eName]);
  ofstream ofs_r2_ev (fname_lst[eGRAPH_ST_R2][eSinc][eName]);
  
  //  FILE* fp_energy = fopen(fname_energy,"w");
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_ST_R1][eSinc][eName]);
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_EV_R1][eSinc][eName]);
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_ST_R2][eSinc][eName]);
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_EV_R2][eSinc][eName]);
  
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
	  sprintf(fname_lst[eZU_E1][eSinc][eName],
		  fname_lst[eZU_E1][eSinc][eTemp], i);
	  sprintf(fname_lst[eZU_E2][eSinc][eName],
		  fname_lst[eZU_E2][eSinc][eTemp], i);
	  
	  ofstream ofs1(fname_lst[eZU_E1][eSinc][eName]);
	  ofstream ofs2(fname_lst[eZU_E2][eSinc][eName]);
	  // ofstream ofs3(fname_lst[eZU_E2][eSinc][eName]);
	  //zu( ofs1 ); all in
	  zu( ofs1, ofs2 );
	  //	zu( ofs1, ofs2, ofs3 ); !not implement
	}
      // output volt ampere field
      // plot("");
    }
}

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

double IDPDIRECT::det(double a11, double a12, double a21, double a22)
{
  return a11*a22 - a12*a21;
}

double IDPDIRECT::En()
{
  double detJ, J11, J12, J21, J22, coef;
  double cth1, cth12, cth2;
  double sth1, sth12, sth2;  
  double PL3, XXa, YYa;
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

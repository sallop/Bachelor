#include <cfloat>
#include "IDP.hpp"

FILE *g_r3;

IDPZCurve::IDPZCurve(){}
IDPZCurve::~IDPZCurve(){}
IDPZCurve::IDPZCurve(string _name,
		     double _t_i, double _t_f,
		     double _dt, double _kiza,
		     int _n_repeat, int _n_divide, int _NN)
  : IDP( _name,	 _t_i,  _t_f, _dt,  _kiza, _n_repeat, _n_divide, _NN)
{//debug();
  printf("%s %d\n", __FUNCTION__, __LINE__);
}

void IDPZCurve::Init(//ScaraRobot _state_i,VRManipulator _v1_i,
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

void IDPZCurve::sinc()
{
  ofstream ofs_energy(     fname_lst[eENERGY][eSinc][eName]);
  ofstream ofs_endeffector(fname_lst[eGRAPH_ENDEFFECTOR][eSinc][eName]);
  ofstream ofs_r1_st (     fname_lst[eGRAPH_ST_R1][eSinc][eName]);
  ofstream ofs_r2_st (     fname_lst[eGRAPH_ST_R2][eSinc][eName]);  
  ofstream ofs_r3_st (     fname_lst[eGRAPH_ST_R3][eSinc][eName]);
  ofstream ofs_r1_ev (     fname_lst[eGRAPH_EV_R1][eSinc][eName]);
  ofstream ofs_r2_ev (     fname_lst[eGRAPH_EV_R2][eSinc][eName]);
  ofstream ofs_r3_ev (     fname_lst[eGRAPH_EV_R3][eSinc][eName]);
  
  vector<Point> trjc;
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_ST_R1][eSinc][eName]);
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_ST_R2][eSinc][eName]);
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_ST_R3][eSinc][eName]);
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_EV_R1][eSinc][eName]);
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_EV_R2][eSinc][eName]);
  printf("line=%d %s\n", __LINE__, fname_lst[eGRAPH_EV_R3][eSinc][eName]);

  g_r3 = fopen("dir-etc/r3zv.dat","w");
  
  double t = 0;
  for( int i = 0; i <= static_cast<int>(Ts/dt); ++i )
    {
      // {{{--------------------------------------------------
      graph_standard(ofs_r1_st, t     ,
       		     r1.th    , r1.thv, r1.tha,
       		     e1.ev    , e1.ia , e1.tau,
       		     e1.energy, e1.ev*e1.ia);
      graph_volt_constituent(ofs_r1_ev, t, e1.ev ,
       			     b1*r1.thv, b2*r1.tha, b3*e1.tau);
      // --------------------------------------------------
      graph_standard(ofs_r2_st, t     ,
       		     r2.th    , r2.thv, r2.tha,
       		     e2.ev    , e2.ia , e2.tau,
       		     e2.energy, e2.ev*e2.ia    );
       
      graph_volt_constituent(ofs_r2_ev, t, e2.ev ,
       			     b1*r2.thv, b2*r2.tha, b3*e2.tau);
      // --------------------------------------------------
      // r3's torque varialbes
//      graph_standard(ofs_r3_st, t,
//       		     r3.z0, r3.Thv(), r3.Tha(),
//       		     e3.ev, e3.ia, e3.tau,
//       		     e3.energy, e3.ev*e3.ia);
//       
//      graph_volt_constituent(ofs_r3_ev, t, e2.ev,
//       			     b1*r3.Thv(), b2*r3.Tha(), b3*e3.tau);
      // }}}--------------------------------------------------
      // r3's endeffector variable
      graph_standard(ofs_r3_st, t    ,
       		     r3.z0    , r3.zv, r3.za ,
       		     e3.ev    , e3.ia, e3.tau,
       		     e3.energy, e3.ev*e3.ia   );
      
      graph_volt_constituent(ofs_r3_ev  , t          , e3.ev    ,
       			     b1*r3.Thv(), b2*r3.Tha(), b3*e3.tau );
      // }}}--------------------------------------------------
      
      // 
      t       = t + dt;
      v1.th   = sc.sin_th1(t);	// virtual machine theta. r angle.
      v1.thv  = sc.sin_th1v(t);
      v1.tha  = sc.sin_th1a(t);
      enemin2 = En();
      {
	double trjc_r = fabs(pnt_f.x - pnt_i.x)/2.0;
	double r2_x2  = trjc_r*trjc_r - v1.x*v1.x;
	fprintf(g_r3,"%8.4lf %8.4lf %8.4lf %8.4lf "
		, t , r3.z0, r3.zv, r3.za);
	fprintf(g_r3,"%8.4lf %8.4lf %8.4lf %8.4lf\n"
		, v1.x, v1.xv, v1.xa, r2_x2);
	//fprintf(fp, "%lf %lf %lf %lf\n", t, v1.x, v1.y);
      }
      ofs_energy << t         <<" "
		 << enemin2   <<" "
		 << e1.energy <<" "
		 << e2.energy <<" "
		 << e3.energy << endl;
      
      graph_endeffector(ofs_endeffector, t, v1, r3);
      
      if (i % 10 == 0)
       	{
	  trjc.push_back( Point(v1.x, v1.y, r3.z0) );

	  
	  sprintf(fname_lst[eZU_E1][eSinc][eName],
		  fname_lst[eZU_E1][eSinc][eTemp], i);
       	  sprintf(fname_lst[eZU_E2][eSinc][eName],
		  fname_lst[eZU_E2][eSinc][eTemp], i);
	  sprintf(fname_lst[eZU_E3][eSinc][eName],
		  fname_lst[eZU_E3][eSinc][eTemp], i);
	  
       	  //for(int i=0; i < eOutputFileNameSize; i++)
	  //  for(int j=0; j < eWhereSize; j++)
	  //    {
	  //	printf("[%2d] ", i);
	  //	printf("eTemp:'%s'\t", fname_lst[i][j][eTemp]);
	  //	printf("eName:'%s'\n", fname_lst[i][j][eName]);		  
	  //    }
	  //assert(false);
	  //printf("i=%d eZU_E1=%d eSinc=%d eName=%d eTemp=%d\n"
	  //, i, eZU_E1, eSinc, eName, eTemp);
 	  //printf(fname_lst[eZU_E1][eSinc][eTemp], i);

	  printf("line=%d %s\n", __LINE__, fname_lst[eZU_E1][eSinc][eName]);
	  
//    	  printf(fname_lst[eZU_E2][eSinc][eTemp], i);
// 	  printf(fname_lst[eZU_E3][eSinc][eTemp], i);
	  
       	  ofstream ofs1(fname_lst[eZU_E1][eSinc][eName]);
       	  ofstream ofs2(fname_lst[eZU_E2][eSinc][eName]);
	  ofstream ofs3(fname_lst[eZU_E3][eSinc][eName]);
       	  //zu( ofs1 ); all in
       	  //zu( ofs1, ofs2 );
	  zu( ofs1, ofs2, ofs3, trjc );// !not implement
       	}
      // output volt ampere field
      // plot("");
    }
  fclose(g_r3);
  // assert(false);
}

void IDPZCurve::sim(int cur, double t)//Magick_numbers
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

inline double sq(double x)
{
  return x*x;
}

inline double cb(double x)
{
  return x*x*x;
}

double IDPZCurve::det(double a11, double a12, double a21, double a22)
{
  return a11*a22 - a12*a21;
}

double IDPZCurve::En()
{
  double detJ, J11, J12, J21, J22, coef;
  double cth1, cth12, cth2;
  double sth1, sth12, sth2;
  static double trjc_r = fabs(pnt_f.x - pnt_i.x)/2.0;
  double PL3, XXa, YYa;
  //  double fxx, fyy;
  
  cth1 = cos(r1.th); cth12 = cos(r1.th + r2.th); cth2 = cos(r2.th);
  sth1 = sin(r1.th); sth12 = sin(r1.th + r2.th); sth2 = sin(r2.th);
  
  double c_vth = cos(v1.th);
  double s_vth = sin(v1.th);
  double t_vth = tan(v1.th);
  //dtan_vth = v1.thv/( cos_vth*cos_vth );
  // ddtan_vth = ( v1.tha*cos_vth + 2*v1.thv*v1.thv*sin_vth )
  // /( cos_vth*cos_vth*cos_vth );
  
  fk1.a = (pnt_f.y - pnt_i.y)/(pnt_f.x - pnt_i.x);// constant
  fk1.b =  pnt_i.y;				  // constant
  fk2.a = tan(v1.th);		// fk2.a = tan(th_.s); variable
  fk2.b = v1.y0 - fk2.a*v1.x0 ; // => check x0         variable
  v1.x  = (fk2.b - fk1.b)/(fk1.a - fk2.a);
  v1.y  =  fk2.a*v1.x + fk2.b;
  
  v1.l  = hypot(v1.x - v1.x0 , v1.y - v1.y0);
  coef  = ( fk1.a*tan(v1.th) + 1.0 )/( fk1.a - tan(v1.th) );
  v1.lv = coef*v1.l*v1.thv;  // v1.lv = coef*v1.l*th_.v;
  v1.xv = v1.lv*cos(v1.th) - v1.l*v1.thv*sin(v1.th);
  v1.yv = v1.lv*sin(v1.th) + v1.l*v1.thv*cos(v1.th);
  
  keirosokudo = sqrt(v1.xv*v1.xv + v1.yv*v1.yv);
  J11 = -r1.l*sth1 - r2.l*sth12;
  J12 = -r2.l*sth12;
  J21 =  r1.l*cth1 + r2.l*cth12;
  J22 =  r2.l*cth12;
  
  detJ   = det(J11, J12, J21, J22);
  r1.thv = det(v1.xv, J12  , v1.yv, J22  )/detJ;
  // => arithmetic exception
  r2.thv = det(J11  , v1.xv, J21  , v1.yv)/detJ;
  //add
  r3.zv = 0.0;			// ???
  
  PL3   = sqrt( v1.x*v1.x + v1.y*v1.y ); // v1.x ? (v1.x - v1.x0)
  r1.th = atan2( v1.y, v1.x ) - acos( 0.5*PL3/r1.l );
  r2.th = 2.0*acos( 0.5*PL3/r1.l );
  
  v1.la  = coef*(v1.l*v1.tha + 2.0*v1.lv*v1.thv) + v1.l*v1.thv*v1.thv;
  v1.xa  = v1.la*cos(v1.th)  - 2.0*v1.lv*v1.thv*sin(v1.th)
    - v1.l*v1.tha*sin(v1.th) - v1.l*v1.thv*v1.thv*cos(v1.th);
  v1.ya  = v1.la*sin(v1.th)  + 2.0*v1.lv*v1.thv*cos(v1.th)
    + v1.l*v1.tha*cos(v1.th) - v1.l*v1.thv*v1.thv*sin(v1.th);
  
  XXa = v1.xa + r1.thv*r1.thv*r1.l*cth1
    + (r1.thv + r2.thv)*(r1.thv + r2.thv)*r2.l*cth12;
  YYa = v1.ya + r1.thv*r1.thv*r1.l*sth1
    + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*sth12;

  // IDPDIRECT modify
  // trjc_r is a radius from init point to end point.
  double r2_x2 = trjc_r*trjc_r - v1.x*v1.x;
  double sqrt_r2_x2 = sqrt(r2_x2);
  //fprintf(stderr, "%g %g\n", r2_x2, DBL_EPSILON);
  //  if (r2_x2 > DBL_EPSILON) {
  if (sqrt_r2_x2 > DBL_EPSILON) {  
    // (dx/dth)
    //double a_tvth = fk1.a - t_vth;
    //double denom[3] = {
    //  v1.y0 - fk1.b - v1.x0*t_vth,
    //  v1.y0 - fk1.b - v1.x0*fk1.a,
    //  2.0*(v1.y0 - fk1.b - v1.x0*fk1.a)*(s_vth*c_vth*(fk1.a - t_vth) + 1.0),
    //};
    //double numer[3] = {
    //  a_tvth,
    //  (c_vth*c_vth)*(a_tvth*a_tvth),
    //  (c_vth*c_vth*c_vth*c_vth)*(a_tvth*a_tvth*a_tvth),
    //};
    //double dzdx[3] = {		// d(z^n)/d(x^n) = z, z', z''
    //  sqrt_r2_x2,		// z = sqrt( r^2 - x^2 );
    //  -v1.x/sqrt_r2_x2,		// dz/dx
    //  -(trjc_r*trjc_r)/(sqrt_r2_x2*sqrt_r2_x2*sqrt_r2_x2), // ddz/dx^2
    //};
    //double dxdth[3]= {
    //  numer[0]/denom[0],
    //  numer[1]/denom[1],
    //  numer[2]/denom[2],
    //};
    //double dthdt[3]= { v1.th, v1.thv, v1.tha};
    
    // (dz/dx)*(dx/dt); dx/dt is v1.xv
    // pderiv 
    r3.z0 = sqrt_r2_x2;
    r3.zv = (-v1.x)*v1.xv/( sqrt_r2_x2 );
    r3.za = (-trjc_r*trjc_r/(sqrt_r2_x2*sqrt_r2_x2*sqrt_r2_x2))*(v1.xv*v1.xv)
      + (-v1.x/sqrt_r2_x2)*(v1.xa);
    
    //r3.za = (-trjc_r*trjc_r)*v1.xv*v1.xv/( r3.z0*r3.z0*r3.z0 )
    // - v1.x*v1.xa/r3.z0;
    
    // deriv 
    //  r3.z0 = sqrt_r2_x2;
    //  r3.zv = -v1.x/sqrt_r2_x2;
    //  r3.za = -trjc_r*trjc_r/r2_x2;

    // chain (dz/dx)*(dx/dth)*(dth/dt)
    // dth is v1.th. dth dt is v1.thv
    // d( (dz/dx)*(dx/dth)*(dth/dt) )/dt
    // r3.z0 = dzdx[0];
    //    r3.zv = dzdx[1]*dxdth[1]*dthdt[1];
    //r3.zv =
    //  -(v1.x/sqrt_r2_x2)
    //  *( (v1.y0 - fk1.b - v1.x0*fk1.a)/(c_vth*c_vth*(fk1.a - t_vth)*(fk1.a - t_vth) ) )
    //  *v1.thv;
    //r3.za
    //  = dzdx[2]*(dxdth[1]*dxdth[1])*(dthdt[1]*dthdt[1])
    //  + dzdx[1]*dxdth[2]*(dthdt[1]*dthdt[1])
    //  + dzdx[1]*dxdth[1]*dthdt[2];
  }
  else {
    r3.z0 = 0.0;
    r3.zv = 0.0;
    r3.za = 0.0;
  }
  
  r1.tha = det(XXa, J12, YYa, J22)/detJ;
  r2.tha = det(J11, XXa, J21, YYa)/detJ;
  
  cth1 = cos(r1.th); cth12 = cos(r1.th + r2.th); cth2 = cos(r2.th);
  sth1 = sin(r1.th); sth12 = sin(r1.th + r2.th); sth2 = sin(r2.th);
  
  e1.tau = (a1+a2+a3+2.0*a4*cth2)*r1.tha + (a2+a3+a4*cth2)*r2.tha 
    - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*sth2;
  e2.tau = (a2 + a3 + a4*cth2)*r1.tha + (a2 + a3)*r2.tha + a4*r1.thv*r1.thv*sth2;
  e3.tau = r3.m*(r3.za + g)*r3.radius;// pulley's radius and torque.

  e1.ev =  b1*r1.thv + b2*r1.tha + b3*e1.tau;
  e1.ia = (e1.ev - Kv*r1.thv)/Ra;
  e2.ev =  b1*r2.thv + b2*r2.tha + b3*e2.tau;
  e2.ia = (e2.ev - Kv*r2.thv)/Ra;
  //  e3.ev = b1*r3.Thv() + b2*r3.Tha() + b3*e3.tau;
  //  e3.ia = (e3.ev - Kv*r3.Thv())/Ra;
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
  
  //if (ev1*ia1 > 0.0){ Ene_1 = Ene_1 + ev1*ia1*dt; }
  //if (ev2*ia2 > 0.0){ Ene_2 = Ene_2 + ev2*ia2*dt; }
#define WATT( e )  ( (e).ev*(e).ia )
  if ( WATT(e1) > 0.0){
    e1.energy = e1.energy + WATT(e1)*dt;
  }
  if ( WATT(e2) > 0.0){ //      Ene_2 = Ene_2 + ev2*ia2*dt;
    e2.energy = e2.energy + WATT(e2)*dt;
  }
  if ( WATT(e3) > 0.0){ //      Ene_2 = Ene_2 + ev2*ia2*dt;
    e3.energy = e3.energy + WATT(e3)*dt;
  }
#undef WATT
  //Ene123 = e1.energy + e2.energy + e3.energy;
  //  return Ene123; // => Nov 23
  //Ene123 = e1.energy + e2.energy;
  return e1.energy + e2.energy + e3.energy;
}

void
IDPZCurve::set_position()
{
  r1.x = r1.x0 + r1.l*cos(r1.th);
  r1.y = r1.y0 + r1.l*sin(r1.th);
  //  r1.z = r1.z0 ;
  r2.x0 = r1.x ; r2.x = r2.x0 + r2.l*cos( r1.th + r2.th );
  r2.y0 = r1.y ; r2.y = r2.y0 + r2.l*sin( r1.th + r2.th );
  r2.z0 = r1.z0;

  r3.x0 = r2.x ;
  r3.y0 = r2.y;
  // r3.z0 = r2.z0;
  r3.z  = r3.z0 + r3.l;
  
  v1.x = v1.x0 + v1.l*cos(v1.th);
  v1.y = v1.y0 + v1.l*sin(v1.th);
  v1.z0 = r3.z0;
}

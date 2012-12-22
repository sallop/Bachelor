#include "IDP.hpp"

void IDPNagato::sinc()
{
  //char fname_energy[64]; sprintf(fname_energy, ftemp_energy, 0);
  //ofstream ofs_energy(fname_energy);
  //ofstream ofs_r1_st(fname_r1_sinc_st), ofs_r1_ev(fname_r1_sinc_ev);
  //ofstream ofs_r2_st(fname_r2_sinc_st), ofs_r2_ev(fname_r2_sinc_ev);
  //ofstream ofs_endeffector(fname_sinc_endeffector);
  ofstream ofs_energy(fname_lst[eENERGY][eSinc][eName]);
  ofstream ofs_r1_st(fname_lst[eGRAPH_ST_R1][eSinc][eName]);
  ofstream ofs_r1_ev(fname_lst[eGRAPH_EV_R1][eSinc][eName]);
  ofstream ofs_r2_st(fname_lst[eGRAPH_ST_R2][eSinc][eName]);
  ofstream ofs_r2_ev(fname_lst[eGRAPH_EV_R2][eSinc][eName]);
  ofstream ofs_endeffector(fname_lst[eGRAPH_ENDEFFECTOR][eSinc][eName]);
  
  //  FILE* fp_energy = fopen(fname_energy,"w");
  //printf("line=%d %s\n", __LINE__, fname_r1_sinc_st);
  //printf("line=%d %s\n", __LINE__, fname_r2_sinc_st);
  //printf("line=%d %s\n", __LINE__, fname_r1_sinc_ev);
  printf("line=%d %s\n", __LINE__
	 , fname_lst[eGRAPH_ST_R1][eSinc][eName]);
  printf("line=%d %s\n", __LINE__
	 , fname_lst[eGRAPH_EV_R1][eSinc][eName]);
  printf("line=%d %s\n", __LINE__
	 , fname_lst[eGRAPH_ST_R2][eSinc][eName]);
  printf("line=%d %s\n", __LINE__
	 , fname_lst[eGRAPH_EV_R2][eSinc][eName]);
  sleep(2);
  
  double t = 0;
  for(nt = 0; nt <= static_cast<int>(Ts/dt); ++nt )
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
      
      graph_endeffector(ofs_endeffector,  t, v1, r3);
      
      t = t + dt;
      v1.th  = sc.sin_th1(t);
      v1.thv = sc.sin_th1v(t);
      v1.tha = sc.sin_th1a(t);

      enemin2 = En();
      // if (nt < nt11){ enemin2 = En1(); }
      // if (nt >= nt11 && nt < nt22){ enemin2 = En2(); }
      // if (nt >= nt22){ enemin2 = En3(); }
      
      ofs_energy << t <<" "<< enemin2 << endl;
      //      fprintf(fp_energy, "%8.4lf %8.4lf\n", t, enemin2 );

      if (nt % 10 == 0)
	{
	  //char fname[4][64];
	  //sprintf(fname[1], ftemp_sinc_e[1], nt);
	  //sprintf(fname[2], ftemp_sinc_e[2], nt);	
	  //sprintf(fname[3], ftemp_sinc_e[3], nt);
	  sprintf(fname_lst[eZU_E1][eSinc][eName],
		  fname_lst[eZU_E1][eSinc][eTemp], nt);
	  sprintf(fname_lst[eZU_E2][eSinc][eName],
		  fname_lst[eZU_E2][eSinc][eTemp], nt);
	  sprintf(fname_lst[eZU_E3][eSinc][eName],
		  fname_lst[eZU_E3][eSinc][eTemp], nt);
	  
	  ofstream ofs1(fname_lst[eZU_E1][eSinc][eName]);
	  ofstream ofs2(fname_lst[eZU_E2][eSinc][eName]);
	  ofstream ofs3(fname_lst[eZU_E3][eSinc][eName]);
	  //zu( ofs1 ); all in
	  //zu( ofs1, ofs2 );
	  zu( ofs1, ofs2, ofs3 );// !not implement
	}
      // output volt ampere field
      // plot("");
    }
  //assert(false);
}

void IDPNagato::set_position()
{
  // using zu()
  r1.x = r1.x0 + r1.l*cos(r1.th);
  r1.y = r1.y0 + r1.l*sin(r1.th);
  r2.x0 = r1.x;
  r2.y0 = r1.y;
  r2.x = r2.x0 + r2.l*cos(r1.th + r2.th);
  r2.y = r2.y0 + r2.l*sin(r1.th + r2.th);
  r3.x0 = r2.x;
  r3.y0 = r2.y;
  r3.z = r3.z0 + r3.l;
  v1.z0 = r3.z0  ;
}

void IDPNagato::sim(int cur, double t)//Magick_numbers
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

IDPNagato::IDPNagato(){};
IDPNagato::~IDPNagato(){};
IDPNagato::IDPNagato(string _name,
			  double _t_i, double _t_f,
			  double _dt, double _kiza,
			  int _n_repeat, int _n_divide, int _NN)
  :IDP( _name, _t_i,  _t_f, _dt,  _kiza, _n_repeat, _n_divide, _NN)
{
  // nt11 = static_cast<int>(0.3125*Ts/dt);
  nt11 = 220;
  nt22 = 440;
  cout << __FILE__ <<" "<<__FUNCTION__ << endl;

//sprintf(fname_r1_sinc_st, "%s/%s_r1_sinc_st.dat"
//	  , dname_dir_etc, name.c_str() );
//sprintf(fname_r2_sinc_st, "%s/%s_r2_sinc_st.dat"
//	  , dname_dir_etc, name.c_str() );
//sprintf(fname_r1_sinc_ev, "%s/%s_r1_sinc_ev.dat"
//	  , dname_dir_etc, name.c_str() );
//sprintf(fname_r2_sinc_ev, "%s/%s_r2_sinc_ev.dat"
//	  , dname_dir_etc, name.c_str() );
//
//sprintf(ftemp_r1_saiteki_st, "%s/%s_r1_saiteki_st_%%05d.dat"
//	  , dname_dir_etc, name.c_str() );
//sprintf(ftemp_r2_saiteki_st, "%s/%s_r2_saiteki_st_%%05d.dat"
//	  , dname_dir_etc, name.c_str() );
//sprintf(ftemp_r1_saiteki_ev, "%s/%s_r1_saiteki_ev_%%05d.dat"
//	  , dname_dir_etc, name.c_str() );
//sprintf(ftemp_r2_saiteki_ev, "%s/%s_r2_saiteki_ev_%%05d.dat"
//	  , dname_dir_etc, name.c_str() );
//
//sprintf(fname_r1_enemin2, "%s/%s_r1_enemin2.dat"
//	  , dname_dir_etc, name.c_str() );
//sprintf(fname_r2_enemin2, "%s/%s_r2_enemin2.dat"
//	  , dname_dir_etc, name.c_str() );
}

void IDPNagato::Init(Model		_state_i,
		     SinCurv		_sc,
		     Linear_curve_t	_fk1,
		     Linear_curve_t	_fk2,
		     Linear_curve_t	_fk3,
		     double		_THs_i, double _THs_f,
		     double		_THv_i, double _THv_f,
		     double		_THa_i, double _THa_f)
		     //Point		pos_i,
		     //Point              pos_f)
{
  IDP::Init(_state_i,_sc,
  	    _THs_i, _THs_f,
  	    _THv_i, _THv_f,
	    _THa_i, _THa_f);
  fk1 = _fk1;
  fk2 = _fk2;
  fk3 = _fk3;
  //ofstream ofs("dir-etc/Init_position.dat");
  //zu(ofs);
  //assert(false);
}

//  double IDPNagato::En(int nt)
double IDPNagato::En()
{
  double ene;
  if (nt < nt11){ ene =  En1(); }
  if (nt >= nt11 && nt < nt22){ ene = En2(); }
  if (nt >= nt22){ ene = En3(); }
  return ene;
}

double IDPNagato::En1()
{
  double Ene12;
  double  fcth1, fcth2, fsth1,fsth2 , fcth12, fsth12;
  double  detJ , XXa, YYa;
  //double  fka1, fkb1,fka2,fkb2, fxv,fyv,PL3; 
  double PL3;
  //fcth1=cos(fth1); fcth12=cos(fth1+fth2); fcth2=cos(fth2); 
  //fsth1=sin(fth1); fsth12=sin(fth1+fth2); fsth2=sin(fth2); 
  fcth1=cos(r1.th); fcth12=cos(r1.th + r2.th); fcth2=cos(r2.th); 
  fsth1=sin(r1.th); fsth12=sin(r1.th + r2.th); fsth2=sin(r2.th); 
  
  // if(nt==1) { frc1= 0.148 ; frc_a=-145.0; }
  
  //Dec28 fxc1=0.0 ; fka1=1.0 ;fkb1=0.05;
  
  // fyc1 = 0.1 ;
	
  //fyc1=fyc1-nt*(0.1-0.06)/(nt11*1.0) ;//*********************2008-2-19
	
  //Dec28 if( nt==1) { fyc1=0.1; }
	
  //Dec28 fyc1= fyc1 - 1.56*(0.04/0.2)*sin(3.14*nt*dt/0.2)*dt ;
  // v1.y0 = v1.y0 - 1.56*(0.04/0.2)*sin(3.14*nt*dt/0.2)*dt;
  // fxc1=0.0; fyc1 =0.06; fka1=1.0 ;fkb1=0.05;
  
  //fka2 = tan(th_1); fkb2=fyc1-fka2*fxc1;
  fk2.a = tan(v1.th);
  fk2.b = v1.y0 - fk2.a*v1.x0;
  
  // printf("%8.4f \n",tan(th_1));
  
  //fxx = (fk2.b - fk1.b)/(fk1.a - fk2.a);
  //fyy = fk1.a*fxx+fk1.b ;
  v1.x = (fk2.b - fk1.b)/(fk1.a - fk2.a);
  v1.y = fk1.a*v1.x + fk1.b ;
  
  //frc_fff=frc1 ;
  
  //frc1=sqrt(pow((fxc1-fxx),2)+pow((fyc1-fyy),2));
  v1.l = hypot( v1.x - v1.x0, v1.y - v1.y0);
  
  //	if(nt==nt11){ printf("nt=%3d  frc1=%8.4f \n",nt,frc1);}
	
  //frc_v= (-frc1*thv_1*(fka1*sin(th_1)+cos(th_1)))/(sin(th_1)-fka1*cos(th_1));
  v1.lv = (-v1.l*v1.thv*(fk1.a*sin(v1.th)+cos(v1.th)))/(sin(v1.th)-fk1.a*cos(v1.th));
  
  //fxv=frc_v*cos(th_1)-frc1*thv_1*sin(th_1);
  //fyv=frc_v*sin(th_1)+frc1*thv_1*cos(th_1);
  v1.xv = v1.lv*cos(v1.th) - v1.l*v1.thv*sin(v1.th);
  v1.yv = v1.lv*sin(v1.th) + v1.l*v1.thv*cos(v1.th);

  // ######
  
  // keirosokudo=sqrt(pow(fxv,2)+pow(fyv,2)) + (1.0/1.4142)*(0.1-0.06)/(nt11*1.0)/dt;

  keirosokudo = hypot(v1.xv, v1.yv);// + 1.56*(0.04/0.22)*(1.0/1.4142)*sin(3.14*nt*dt/0.22);
  
  // keirosokudo=-frc1*thv_1;

  //detJ=L1*L2*fsth2;
  detJ = r1.l*r2.l*fsth2;
  
  //fthv1=(fxv*L2*fcth12+fyv*L2*fsth12)/detJ ;
  //fthv2=(-fyv*(L1*fsth1+L2*fsth12)-fxv*(L1*fcth1+L2*fcth12))/detJ ;
  r1.thv = (v1.xv*r2.l*fcth12 + v1.yv*r2.l*fsth12)/detJ;
  r2.thv = (-r2.yv*(r2.l*fsth1 + r2.l*fsth12) - v1.xv*(r1.l*fcth1 + r2.l*fcth12))/detJ;
  
  //PL3=sqrt(fxx*fxx + fyy*fyy);
  PL3 = hypot( v1.x - r1.x0, v1.y - r1.y0);	// ?
  
  //fth1=atan2(fyy,fxx)-acos(0.5*PL3/L1);
  //fth2=2*acos(0.5*PL3/L1) ;
  r1.th = atan2(v1.y, v1.x) - acos(0.5*PL3/r1.l);
  r2.th = 2.0*acos(0.5*PL3/r1.l);

  
  
  //frc_a=((-2*frc_v*thv_1-frc1*tha_1)*(fka1*sin(th_1)+cos(th_1))-frc1*thv_1*thv_1*(fka1*cos(th_1)-sin(th_1)))/(sin(th_1)-fka1*cos(th_1));
  
  //frc_a =frc1*thv_1*thv_1+(2*frc_v*thv_1+tha_1*frc1)*(fka1*sin(th_1)+cos(th_1))/(fka1*cos(th_1)-sin(th_1));
  v1.la = v1.l*v1.thv*v1.thv
    + (2.0*v1.lv*v1.thv + v1.tha*v1.l)*(fk1.a*sin(v1.th)+cos(v1.th))/(fk1.a*cos(v1.th)-sin(v1.th));
  
  //fxa=frc_a*cos(th_1)-2*frc_v*thv_1*sin(th_1)-frc1*tha_1*sin(th_1)-frc1*pow(thv_1,2)*cos(th_1);
  //fya=frc_a*sin(th_1)+2*frc_v*thv_1*cos(th_1)+frc1*tha_1*cos(th_1)-frc1*pow(thv_1,2)*sin(th_1);
  v1.xa = v1.la*cos(v1.th) - 2.0*v1.lv*v1.thv*sin(v1.th)
    - v1.l*v1.tha*sin(v1.th) - v1.l*v1.thv*v1.thv*cos(v1.th);
  v1.ya = v1.la*sin(v1.th) + 2*v1.lv*v1.thv*cos(v1.th) + v1.l*v1.tha*cos(v1.th)
    - v1.l*v1.thv*v1.thv*sin(v1.th);
  
  //XXa = fxa + pow(fthv1,2)*L1*fcth1 + pow((fthv1+fthv2),2)*L2*fcth12;
  //YYa = fya + pow(fthv1,2)*L1*fsth1 + pow((fthv1+fthv2),2)*L2*fsth12;
  XXa = v1.xa + r1.thv*r1.thv*r1.l*fcth1
    + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*fcth12;
  YYa = v1.ya + r1.thv*r1.thv*r1.l*fsth1
    + (r1.thv + r2.thv)*(r1.thv + r2.thv)*r2.l*fsth12;
  
  //ftha1=(L2*fcth12*XXa  +  L2*fsth12*YYa)/detJ ;
  //ftha2=((-L1*fcth1-L2*fcth12)*XXa + (-L1*fsth1-L2*fsth12)*YYa)/detJ ;
  r1.tha = (r2.l*fcth12*XXa + r2.l*fsth12*YYa)/detJ;
  r2.tha = ((-r1.l*fcth1 - r2.l*fcth12)*XXa + (-r1.l*fsth1 - r2.l*fsth12)*YYa)/detJ;
  
  //fcth1=cos(fth1); fcth12=cos(fth1+fth2);  fcth2=cos(fth2); 
  //fsth1=sin(fth1); fsth12=sin(fth1+fth2);  fsth2=sin(fth2);
  fcth1 = cos(r1.th); fcth12 = cos(r1.th + r2.th);  fcth2 = cos(r2.th); 
  fsth1 = sin(r1.th); fsth12 = sin(r1.th + r2.th);  fsth2 = sin(r2.th); 
  //---------- 2 Link System -------------------------------------
  //tau_1=(a1+a2+2*a3*fcth2)*ftha1  + (a2+a3*fcth2)*ftha2 ;
  //tau_1=tau_1 - a3*(2*fthv1*fthv2 + pow(fthv2,2))*fsth2 + a4*fcth1 + a5*fcth12;
  e1.tau = (a1 + a2 + a3 + 2.0*a4*fcth2)*r1.tha + (a2 + a3 + a4*fcth2)*r2.tha
    - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*fsth2;// + a4*fcth1 + a5*fcth12;
  
  e2.tau = (a2 + a3 + a4*fcth2)*r1.tha
    + (a2 + a3)*r2.tha
    + a4*r1.thv*r1.thv*fsth2;// + a5*fcth12;
  
  //ev1 = b1*fthv1  +b2*ftha1  +b3*tau_1 ; ia1=(ev1-Kv*fthv1)/Ra ;
  //ev2 = b1*fthv2  +b2*ftha2  +b3*tau_2 ; ia2=(ev2-Kv*fthv2)/Ra ;
  e1.ev = b1*r1.thv + b2*r1.tha + b3*e1.tau;
  e1.ia = (e1.ev - Kv*r1.thv)/Ra;
  e2.ev = b1*r2.thv + b2*r2.tha + b3*e2.tau;
  e2.ia = (e2.ev - Kv*r2.thv)/Ra;

  //if( ev1*ia1 > 0.0 ){ Ene_1 = Ene_1 + ev1*ia1*dt ; }
  //if( ev2*ia2 > 0.0 ){ Ene_2 = Ene_2 + ev2*ia2*dt ; }
  if( e1.ev*e1.ia > 0.0 ){ e1.energy = e1.energy + e1.ev*e1.ia*dt; }
  if( e2.ev*e2.ia > 0.0 ){ e2.energy = e2.energy + e2.ev*e2.ia*dt; }
  
  //Ene12  = Ene_1 + Ene_2 ;
  Ene12  = e1.energy + e2.energy;
  // Ene12  = e1.energy + e2.energy + e3.energy;
  return e1.energy + e2.energy;
  //return e1.energy + e2.energy + e3.energy;
}

//************************************************************
double IDPNagato::En2()
{
  static const double fxc1 = 0.0, fyc1 = 0.06;
  static const double frc1 = 0.01/sqrt(2.0);
  //static const double frc1 = fyc1/sqrt(2.0);
  double Ene12;
  double  fcth1, fcth2, fsth1,fsth2 , fcth12, fsth12;
  double  detJ , XXa, YYa;
  //double  fka1, fkb1,fka2,fkb2, fxv,fyv;
  double PL3; 
  
  
  //fcth1=cos(fth1); fcth12=cos(fth1+fth2); fcth2=cos(fth2); 
  //fsth1=sin(fth1); fsth12=sin(fth1+fth2); fsth2=sin(fth2); 

  fcth1=cos(r1.th); fcth12=cos(r1.th+r2.th); fcth2=cos(r2.th); 
  fsth1=sin(r1.th); fsth12=sin(r1.th+r2.th); fsth2=sin(r2.th); 

  //if(nt==nt11+1) { frc1= 0.148 ; frc_a=-145.0; } 

  //fxc1=0.0; fyc1 =0.1; //fka1=0.0 ;fkb1=0.1;

  //Dec 28 fxc1=0.0; fyc1 =0.06;//********************************
  // v1.x0 = 0.0;
  // v1.y0 = 0.06;
  //frc1=fyc1/sqrt(2.0) ;
	
  //  frc1=0.0078 ; frc1=0.0365 ;
  //Dec28 v1.l = v1.y0/sqrt(2.0);
	
  //Dec 28 frc1=0.04246 ;//**********************************
	
  //Dec 28 frc1=0.01/sqrt(2.0);

  
  //frc1=0.0354;
  //  frc1=0.038 ; frc1=0.0389 ;
  //frc1= 0.148 ;
  
  // fxx = fxc1 + frc1*cos(th_1); 
  // fyy = fyc1 + frc1*sin(th_1);
  // v1.x = v1.x0 + v1.l*cos(v1.th);
  // v1.y = v1.y0 + v1.l*sin(v1.th);
  v1.x = fxc1 + frc1*cos(v1.th);
  v1.y = fyc1 + frc1*sin(v1.th);
  
  //frc_v= (-frc1*thv_1*cos(th_1))/(sin(th_1));
  //v1.lv= (-v1.l*v1.thv*cos(v1.th))/(sin(v1.th));
  v1.lv= (-frc1*v1.thv*cos(v1.th))/(sin(v1.th));
  
  //fxv=-thv_1*frc1*sin(th_1); 
  //fyv= thv_1*frc1*cos(th_1);
  //v1.xv = -v1.thv*v1.l*sin(v1.th); 
  //v1.yv =  v1.thv*v1.l*cos(v1.th);
  v1.xv = -v1.thv*frc1*sin(v1.th); 
  v1.yv =  v1.thv*frc1*cos(v1.th);
  
  //keirosokudo=sqrt(pow(fxv,2)+pow(fyv,2));
  //keirosokudo=-frc1*thv_1;
  //keirosokudo = -v1.l*v1.thv;
  keirosokudo = -frc1*v1.thv;
	
  //detJ=L1*L2*fsth2;
  detJ=r1.l*r2.l*fsth2;
  
  //fthv1=(fxv*L2*fcth12 + fyv*L2*fsth12)/detJ ;
  //fthv2=(-fyv*(L1*fsth1+L2*fsth12)-fxv*(L1*fcth1+L2*fcth12))/detJ ;
  //r1.thv = (v1.xv*r2.l*fcth12+v1.yv*r2.l*fsth12)/detJ ;
  //r2.thv =(-v1.yv*(r1.l*fsth1+r2.l*fsth12) - v1.xv*(r1.l*fcth1+r2.l*fcth12))/detJ ;
  r1.thv = (v1.xv*r2.l*fcth12 + v1.yv*r2.l*fsth12)/detJ ;
  r2.thv = (-v1.yv*(r1.l*fsth1+r2.l*fsth12) - v1.xv*(r1.l*fcth1+r2.l*fcth12))/detJ ;
  
  //PL3=sqrt(fxx*fxx + fyy*fyy);
  PL3 = hypot(v1.x - r1.x0, v1.y - r1.y0);	// ?
  
  //fth1=atan2(fyy,fxx)-acos(0.5*PL3/L1);
  //fth2=2*acos(0.5*PL3/L1) ;
  r1.th = atan2(v1.y - r1.y0, v1.x - r1.x0) - acos(0.5*PL3/r1.l);
  r2.th = 2.0*acos(0.5*PL3/r1.l);
  
  //frc_a=((-2*frc_v*thv_1-frc1*tha_1)*(cos(th_1))-frc1*thv_1*thv_1*(-sin(th_1)))/(sin(th_1));
  //fxa =-tha_1*frc1*sin(th_1)-thv_1*thv_1*frc1*cos(th_1);
  //fya = tha_1*frc1*cos(th_1)-thv_1*thv_1*frc1*sin(th_1);
  
  v1.la = ((-2.0*v1.lv*v1.thv - frc1*v1.tha)*(cos(v1.th))
	   - frc1*v1.thv*v1.thv*(-sin(v1.th)))/sin(v1.th);
  v1.xa = -v1.tha*frc1*sin(v1.th) - v1.thv*v1.thv*frc1*cos(v1.th);
  v1.ya =  v1.tha*frc1*cos(v1.th) - v1.thv*v1.thv*frc1*sin(v1.th);
  
  //XXa = fxa + pow(fthv1,2)*L1*fcth1 + pow((fthv1+fthv2),2)*L2*fcth12;
  //YYa = fya + pow(fthv1,2)*L1*fsth1 + pow((fthv1+fthv2),2)*L2*fsth12;
  XXa = v1.xa + r1.thv*r1.thv*r1.l*fcth1
      + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*fcth12;
  YYa = v1.ya + r1.thv*r1.thv*r1.l*fsth1 + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*fsth12;
  
  //ftha1=(L2*fcth12*XXa  +  L2*fsth12*YYa)/detJ ;
  //ftha2=((-L1*fcth1-L2*fcth12)*XXa + (-L1*fsth1-L2*fsth12)*YYa)/detJ;
  r1.tha = (r2.l*fcth12*XXa  +  r2.l*fsth12*YYa)/detJ ;
  r2.tha = ((-r1.l*fcth1 - r2.l*fcth12)*XXa + (-r1.l*fsth1 - r2.l*fsth12)*YYa)/detJ ;

  // using zu()
  r1.x = r1.x0 + r1.l*cos(r1.th);
  r1.y = r1.y0 + r1.l*sin(r1.th);
  r2.x0 = r1.x;
  r2.y0 = r1.y;
  r2.x = r2.x0 + r2.l*cos(r1.th + r2.th);
  r2.y = r2.y0 + r2.l*sin(r1.th + r2.th);
  
  //fcth1=cos(fth1); fcth12=cos(fth1+fth2);  fcth2=cos(fth2); 
  //fsth1=sin(fth1); fsth12=sin(fth1+fth2);  fsth2=sin(fth2);
  fcth1=cos(r1.th); fcth12=cos(r1.th+r2.th);  fcth2=cos(r2.th); 
  fsth1=sin(r1.th); fsth12=sin(r1.th+r2.th);  fsth2=sin(r2.th); 
  //---------- 2 Link System -------------------------------------
  //tau_1=(a1+a2+2*a3*fcth2)*ftha1  + (a2+a3*fcth2)*ftha2 ;
  //tau_1=tau_1 - a3*(2*fthv1*fthv2 + pow(fthv2,2))*fsth2 + a4*fcth1 + a5*fcth12;
  e1.tau = (a1 + a2 + a3 + 2.0*a4*fcth2)*r1.tha + (a2 + a3 + a4*fcth2)*r2.tha
    - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*fsth2;// + a4*fcth1 + a5*fcth12;
  
  //tau_2=(a2+a3*fcth2)*ftha1 + a2*ftha2 + a3*pow(fthv1,2)*fsth2 + a5*fcth12;
  e2.tau = (a2 + a3 + a4*fcth2)*r1.tha
    + (a2 + a3)*r2.tha
    + a4*r1.thv*r1.thv*fsth2;// + a5*fcth12;

  //e3 = ;
  
  //ev1 = b1*fthv1  +b2*ftha1  +b3*tau_1 ; ia1=(ev1-Kv*fthv1)/Ra ;
  //ev2 = b1*fthv2  +b2*ftha2  +b3*tau_2 ; ia2=(ev2-Kv*fthv2)/Ra ;
  e1.ev = b1*r1.thv + b2*r2.tha + b3*e1.tau;
  e1.ia =(e1.ev - Kv*r1.thv)/Ra ;
  
  e2.ev = b1*r2.thv + b2*r2.tha + b3*e2.tau;
  e2.ia = (e2.ev - Kv*r2.thv)/Ra ; 

  //e3.ev = ;
  //e3.ia = ;
  
  //if( ev1*ia1 > 0.0 ){ Ene_1 = Ene_1 + ev1*ia1*dt ; }
  //if( ev2*ia2 > 0.0 ){ Ene_2 = Ene_2 + ev2*ia2*dt ; }
  if( e1.ev*e1.ia > 0.0 ){ e1.energy = e1.energy + e1.ev*e1.ia*dt ; }
  if( e2.ev*e2.ia > 0.0 ){ e2.energy = e2.energy + e2.ev*e2.ia*dt ; }
  //if( e3.ev*e3.ia > 0.0 ){ e3.energy = e3.energy + e3.ev*e3.ia*dt ; }
  
  //Ene12  = Ene_1 + Ene_2 ;
  Ene12  = e1.energy + e2.energy;
  //Ene12  = e1.energy + e2.energy + e3.energy;
  return e1.energy + e2.energy;
  //return e1.energy + e2.energy + e3.energy;
}

//******************************************************************

double IDPNagato::En3()
{
  double Ene12;
  double  fcth1, fcth2, fsth1,fsth2 , fcth12, fsth12;
  double  detJ , XXa, YYa;
  
  //double  fka1, fkb1,fka2,fkb2, fxv,fyv,PL3; 
  double PL3;
  //fcth1=cos(fth1); fcth12=cos(fth1+fth2); fcth2=cos(fth2); 
  //fsth1=sin(fth1); fsth12=sin(fth1+fth2); fsth2=sin(fth2);
  fcth1=cos(r1.th); fcth12=cos(r1.th + r2.th); fcth2=cos(r2.th); 
  fsth1=sin(r1.th); fsth12=sin(r1.th + r2.th); fsth2=sin(r2.th); 
  
  // if(nt==nt22+1) { frc1= 0.148 ; frc_a=-145.0; } 
  
  //fxc1=0.0; fyc1 =0.1; fka1=-1.0 ;fkb1=0.05;
  
  //Dec28 fxc1=0.0;
  //Dec28 fka1=-1.0;
  //Dec28 fkb1=0.05; //***********
  fk3.a;
  fk3.b;
  
  //fyc1 = 0.06 ;
	
  //if(nt==nt22+1) { fyc1 = 0.06 ; }
  //Dec28 if(nt == nt22+1) { v1.y0 = 0.06 ; }
  
  //fyc1=fyc1 + (nt-nt22)*(0.1-0.06)/((640-nt22)*1.0) ;
	
  //fyc1= fyc1 + 1.56*(0.04/0.2)*sin(3.14*(nt-nt22)*dt/0.2)*dt ;
  //  v1.y0 = v1.y0 + 1.56*(0.04/0.2)*sin(3.14*(nt-nt22)*dt/0.2)*dt ;
  
  //  fka2=tan(th_1); fkb2=fyc1-fka2*fxc1;
  fk2.a = tan(v1.th);
  fk2.b = v1.y0 - fk2.a*v1.x0;
  
  // printf("%8.4f \n",tan(th_1));
  
  //fxx=(fkb2-fkb1)/(fka1-fka2);
  //fyy= fka1*fxx+fkb1 ;
  v1.x = (fk2.b - fk3.b)/(fk3.a - fk2.a);
  v1.y =  fk3.a*v1.x + fk3.b;//?
  
  //frc_fff=frc1 ;
  
  //frc1=sqrt(pow((fxc1-fxx),2)+pow((fyc1-fyy),2));
  v1.l = hypot( v1.x - v1.x0, v1.y - v1.y0);
  
  //frc_v= (-frc1*thv_1*(fka1*sin(th_1)+cos(th_1)))/(sin(th_1)-fka1*cos(th_1));
  v1.lv = (-v1.l*v1.thv*(fk3.a*sin(v1.th)+cos(v1.th)))/(sin(v1.th)-fk3.a*cos(v1.th));
  
  //fxv=frc_v*cos(th_1) - frc1*thv_1*sin(th_1);
  //fyv=frc_v*sin(th_1) + frc1*thv_1*cos(th_1);
  v1.xv = v1.lv*cos(v1.th) - v1.l*v1.thv*sin(v1.th);
  v1.yv = v1.lv*sin(v1.th) + v1.l*v1.thv*cos(v1.th);
  
  //keirosokudo=sqrt(pow(fxv,2)+pow(fyv,2)) +(1.0/1.4142)*(0.1-0.06)/((640-nt22)*1.0)/dt;
  
  //keirosokudo=sqrt(pow(fxv,2)+pow(fyv,2)) + 1.56*(1.0/1.4142)*(0.04/0.2)*sin(3.14*(nt-nt22)*dt/0.2) ;
  keirosokudo = hypot(v1.xv, v1.yv);// + 1.56*(1.0/1.4142)*(0.04/0.2)*sin(3.14*(nt-nt22)*dt/0.2);
  //  keirosokudo=-frc1*thv_1;
  
  //detJ=L1*L2*fsth2;
  detJ = r1.l*r2.l*fsth2;
  
  //fthv1=(fxv*L2*fcth12+fyv*L2*fsth12)/detJ ;
  //fthv2=(-fyv*(L1*fsth1+L2*fsth12)-fxv*(L1*fcth1+L2*fcth12))/detJ ;
  r1.thv = ( v1.xv*r2.l*fcth12 + v1.yv*r2.l*fsth12)/detJ ;
  r2.thv = (-v1.yv*(r1.l*fsth1+r2.l*fsth12) - v1.xv*(r1.l*fcth1 + r2.l*fcth12))/detJ ;
  
  //PL3=sqrt(fxx*fxx + fyy*fyy);
  PL3 = hypot(v1.x - r1.x0, v1.y - r1.y0);
  
  //fth1=atan2(fyy,fxx)-acos(0.5*PL3/L1);
  //fth2=2*acos(0.5*PL3/L1) ;
  r1.th = atan2(v1.y, v1.x) - acos(0.5*PL3/r1.l);
  r2.th = 2.0*acos(0.5*PL3/r1.l);
  
  //frc_a=((-2*frc_v*thv_1-frc1*tha_1)*(fka1*sin(th_1)+cos(th_1))
  //-frc1*thv_1*thv_1*(fka1*cos(th_1)-sin(th_1)))/(sin(th_1)-fka1*cos(th_1));
  
  //frc_a =frc1*thv_1*thv_1
  //+(2*frc_v*thv_1+tha_1*frc1)*(fka1*sin(th_1)+cos(th_1))/(fka1*cos(th_1)-sin(th_1));
  v1.la = v1.l*v1.thv*v1.thv
    +(2.0*v1.lv*v1.thv+v1.tha*v1.l)*(fk3.a*sin(v1.th)+cos(v1.th))/(fk3.a*cos(v1.th)-sin(v1.th));
  
  //fxa=frc_a*cos(th_1)-2*frc_v*thv_1*sin(th_1)-frc1*tha_1*sin(th_1)-frc1*pow(thv_1,2)*cos(th_1);
  //fya=frc_a*sin(th_1)+2*frc_v*thv_1*cos(th_1)+frc1*tha_1*cos(th_1)-frc1*pow(thv_1,2)*sin(th_1);
  v1.xa = v1.la*cos(v1.th) -2.0*v1.lv*v1.thv*sin(v1.th)
    - v1.l*v1.tha*sin(v1.th) - v1.l*v1.thv*v1.thv*cos(v1.th);
  v1.ya = v1.la*sin(v1.th)+2.0*v1.lv*v1.thv*cos(v1.th)+v1.l*v1.tha*cos(v1.th)
    -v1.l*v1.thv*v1.thv*sin(v1.th);
  
  //XXa = fxa + pow(fthv1,2)*L1*fcth1 + pow((fthv1+fthv2),2)*L2*fcth12;
  //YYa = fya + pow(fthv1,2)*L1*fsth1 + pow((fthv1+fthv2),2)*L2*fsth12;
  XXa = v1.xa + v1.thv*v1.thv*r1.l*fcth1 + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*fcth12;
  YYa = v1.ya + v1.thv*v1.thv*r1.l*fsth1 + (r1.thv+r2.thv)*(r1.thv+r2.thv)*r2.l*fsth12;
  
  //ftha1=(L2*fcth12*XXa  +  L2*fsth12*YYa)/detJ ;
  //ftha2=((-L1*fcth1-L2*fcth12)*XXa + (-L1*fsth1-L2*fsth12)*YYa)/detJ ;
  r1.tha = (r2.l*fcth12*XXa  +  r2.l*fsth12*YYa)/detJ ;
  r2.tha = ( (-r1.l*fcth1 - r2.l*fcth12)*XXa + (-r1.l*fsth1 - r2.l*fsth12)*YYa)/detJ ;

  //fcth1=cos(fth1); fcth12=cos(fth1+fth2);  fcth2=cos(fth2); 
  //fsth1=sin(fth1); fsth12=sin(fth1+fth2);  fsth2=sin(fth2);
  fcth1=cos(r1.th); fcth12=cos(r1.th + r2.th);  fcth2=cos(r2.th); 
  fsth1=sin(r1.th); fsth12=sin(r1.th + r2.th);  fsth2=sin(r2.th); 

  //---------- 2 Link System -------------------------------------

  //tau_1=(a1+a2+2*a3*fcth2)*ftha1  + (a2+a3*fcth2)*ftha2 ;
  //tau_1=tau_1 - a3*(2*fthv1*fthv2 + pow(fthv2,2))*fsth2 + a4*fcth1 + a5*fcth12;
  e1.tau = (a1 + a2 + a3 + 2.0*a4*fcth2)*r1.tha + (a2 + a3 + a4*fcth2)*r2.tha
    - a4*(2.0*r1.thv*r2.thv + r2.thv*r2.thv)*fsth2;// + a4*fcth1 + a5*fcth12;
  
  //tau_2=(a2+a3*fcth2)*ftha1 + a2*ftha2 + a3*pow(fthv1,2)*fsth2 + a5*fcth12;
  e2.tau = (a2 + a3 + a4*fcth2)*r1.tha
    + (a2 + a3)*r2.tha
    + a4*r1.thv*r1.thv*fsth2;// + a5*fcth12;  

  //ev1 = b1*fthv1  +b2*ftha1  +b3*tau_1 ; ia1=(ev1-Kv*fthv1)/Ra ;
  //ev2 = b1*fthv2  +b2*ftha2  +b3*tau_2 ; ia2=(ev2-Kv*fthv2)/Ra ;
  e1.ev = b1*r1.thv + b2*r1.tha + b3*e1.tau;
  e1.ia = (e1.ev - Kv*r1.thv)/Ra ;
  e2.ev = b1*r2.thv + b2*r2.tha + b3*e2.tau;
  e2.ia = (e2.ev - Kv*r2.thv)/Ra ; 

  //if( ev1*ia1 > 0.0 ){ Ene_1 = Ene_1 + ev1*ia1*dt ; }
  //if( ev2*ia2 > 0.0 ){ Ene_2 = Ene_2 + ev2*ia2*dt ; }
  if( e1.ev*e1.ia > 0.0 ){ e1.energy = e1.energy + e1.ev*e1.ia*dt; }
  if( e2.ev*e2.ia > 0.0 ){ e2.energy = e2.energy + e2.ev*e2.ia*dt; }
  //if( e2.ev*e2.ia > 0.0 ){ e2.energy = e2.energy + e2.ev*e2.ia*dt; }  
  //Ene12  = Ene_1 + Ene_2 ;
  Ene12  = e1.energy + e2.energy;
  // Ene12  = e1.energy + e2.energy + e3.energy;
  return e1.energy + e2.energy;
  //  return e1.energy + e2.energy + e3.energy;
}

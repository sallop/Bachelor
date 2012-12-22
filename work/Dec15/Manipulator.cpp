#include <iostream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include "Manipulator.hpp"

Manipulator::Manipulator()
{
  m  = 0.0; l  = 0.0; r  = 0.0;
  x0 = 0.0; y0 = 0.0; z0 = 0.0;
//  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<<__LINE__<<std::endl;
//  std::cout << "m=" << m << " l=" << l << " r=" << r << std::endl;
//  std::cout << "x0=" << x0 << " y0=" << y0 << " z0=" << z0 << std::endl;
}

Manipulator::Manipulator(double _m, double _l,
			 double _x0, double _y0, double _z0)
{
  m  = _m ;
  l   = _l;
  r   = _l / 2.0;
  x0 = _x0;
  y0 = _y0;
  z0 = _z0;
//   std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<<__LINE__<<std::endl;
//   std::cout << "m=" << m << " l=" << l << " r=" << r << std::endl;
//  std::cout << "x0=" << x0 << " y0=" << y0 << " z0=" << z0 << std::endl;
}

Manipulator::~Manipulator()
{
  //  std::cout << __FILE__<<"_"<<__FUNCTION__<<"_"<<__LINE__<<std::endl;
}

void Manipulator::print_value(std::ostream& os)
{
  os  << m <<" "<< l <<" "<< r <<" "
      << x0<<" "<< y0<<" "<< z0<< std::endl;
}

void Manipulator::setXYZ0(double _x0, double _y0, double _z0)
{
  x0 = _x0; y0 = _y0; z0 = _z0;
}

void Manipulator::setXYZ(double _x, double _y, double _z)
{
  //  std::cout <<__FILE__<<"_"<<__FUNCTION__<<"_"<<__LINE__<<std::endl;
}


VRManipulator::VRManipulator()
  : Manipulator()
{
  //  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__<< std::endl;
  x  = 0.0; xv  = 0.0; xa  = 0.0;
  y  = 0.0; yv  = 0.0; ya  = 0.0;
  th = 0.0; thv = 0.0; tha = 0.0;
  l  = 0.0; lv  = 0.0; la  = 0.0;
	
  //  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__ << std::endl;
}

VRManipulator::VRManipulator(double _m, double _l,
			     double _x0, double _y0, double _z0,
			     double _rth)
  : Manipulator( _m,  _l,  _x0,  _y0,  _z0)//	: Manipulator()
{
  x  = 0.0; xv  = 0.0; xa  = 0.0;
  y  = 0.0; yv  = 0.0; ya  = 0.0;
  th =_rth; thv = 0.0; tha = 0.0;
            lv  = 0.0; la  = 0.0;
	
	    //  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__ << std::endl;
}

VRManipulator::VRManipulator(const VRManipulator &vrm)
{
  *this = vrm;
}

void VRManipulator::setXYZ0(double _x0, double _y0, double _z0)
{
  x0 = _x0; y0 = _y0; z0 = _z0;
}

void VRManipulator::setXYZ(double _x, double _y, double _z)
{
  x = _x; y = _y; //z = _z;
}

// Dec 2 VRManipulator& VRManipulator::operator=(const VRManipulator& vrm)
// Dec 2 {
// Dec 2   this->m  = vrm.m ;
// Dec 2   this->l  = vrm.l ;
// Dec 2   this->r  = vrm.r ;
// Dec 2   this->x0 = vrm.x0;
// Dec 2   this->y0 = vrm.y0;
// Dec 2   this->z0 = vrm.z0;
// Dec 2 	
// Dec 2   this->x  = vrm.x ; this->xv  = vrm.xv ; this->xa  = vrm.xa ;
// Dec 2   this->y  = vrm.y ; this->yv  = vrm.yv ; this->ya  = vrm.ya ;
// Dec 2   this->th = vrm.th; this->thv = vrm.thv; this->tha = vrm.tha;
// Dec 2   //this->l  = vrm.l ; this->lv  = vrm.lv ; this->la  = vrm.la ;
// Dec 2                      this->lv  = vrm.lv ; this->la  = vrm.la ;  
// Dec 2   //	this->z  = vrm.z ;
// Dec 2   return *this;
// Dec 2 }

VRManipulator::~VRManipulator()
{
  //  std::cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<<std::endl;
}
				 

Revolute::Revolute() : Manipulator()
{
  //  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__ << std::endl;
}

Revolute::Revolute(double _m, double _l,
		   double _x0, double _y0, double _z0, double _rth)
  : Manipulator( _m,  _l,  _x0,  _y0,  _z0)
{
  Ig = (m*l*l) / 12.0;
  th = _rth; thv = 0.0;	tha = 0.0;
  x  =  0.0;  xv = 0.0;   xa = 0.0;
  y  =  0.0;  yv = 0.0;   ya = 0.0;
  //  std::cout << "Ig=" << Ig << std::endl;
  //  std::cout << "th=" << th << std::endl;
}

Revolute::Revolute(const Revolute &revo)
{
  //  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__ << std::endl;
  //  std::cout << "Copy constructor called." << std::endl;
  *this = revo;
}

// Dec 2 Revolute& Revolute::operator=(const Revolute& revo)
// Dec 2 {
// Dec 2   std::cout << __FILE__<<"_"<< __FUNCTION__ <<"_"<< __LINE__ <<std::endl;
// Dec 2   // super class
// Dec 2   this->m = revo.m;
// Dec 2   this->l = revo.l;
// Dec 2   this->r = revo.r;
// Dec 2   this->x0 = revo.x0;
// Dec 2   this->y0 = revo.y0;		
// Dec 2   this->z0 = revo.z0;	
// Dec 2   //
// Dec 2   this->Ig = revo.Ig;
// Dec 2 	
// Dec 2   this->x  = revo.x ;
// Dec 2   this->xv = revo.xv;
// Dec 2   this->xa = revo.xa;
// Dec 2 	
// Dec 2   this->y  = revo.y ;
// Dec 2   this->yv = revo.yv;
// Dec 2   this->ya = revo.ya;
// Dec 2 
// Dec 2   this->th  = revo.th ;
// Dec 2   this->thv = revo.thv;
// Dec 2   this->tha = revo.tha;
// Dec 2 	
// Dec 2   return *this;
// Dec 2 }

Revolute::~Revolute()
{
  //  std::cout << __FILE__ <<"_"<< __FUNCTION__ << "_" << __LINE__ << std::endl;
}

void Revolute::setXYZ0(double _x0, double _y0, double _z0)
{
  x0 = _x0; y0 = _y0; z0 = _z0;
}

void Revolute::setXYZ(double _x, double _y, double _z)
{
  x = _x; y = _y;// z = _z;
}

//
// Th s, v, a is modify,
// but not modify xv, yv, ya
void Revolute::setTHsva(double s, double v, double a)
{
  th = s; thv = v; tha = a; 
}

void Revolute::setXsva(double s, double v, double a)
{
  x = s; xv = v; xa = a; 
}

void Revolute::setYsva(double s, double v, double a)
{
  y = s; yv = v; ya = a; 
}

void Revolute::setEndEffector(double absTh)
{
  x = x0 + l*cos(absTh);
  y = y0 + l*sin(absTh);
	
//   std::cout << __FILE__     <<"_"
// 	    << __FUNCTION__ <<"_"
// 	    << __LINE__     << std::endl;
//   std::cout << "x= " << x << std::endl;
//   std::cout << "y= " << y << std::endl;
}

Prismatic::Prismatic() : Manipulator()
{
//   std::cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<< std::endl;
}

Prismatic::Prismatic(double _m, double _l,
		     double _x0, double _y0, double _z0, double _radius)
  : Manipulator( _m,  _l,  _x0,  _y0,  _z0)
{
  radius = _radius;
  Iz = 0.0;
  z = z0 - l; zv = 0.0; za = 0.0;
		
//   std::cout	<< "radius=" << radius << " Iz="<< Iz
// 		<< " z=" << z << " zv="<< zv << " za="<< za << std::endl;
}

Prismatic::Prismatic(const Prismatic &pris)
{
  //  std::cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<< std::endl;
  *this = pris;
}

Prismatic::~Prismatic()
{
  //  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_" << __LINE__ << std::endl;
}

double Prismatic::Thv(){ return zv/radius;}
double Prismatic::Tha(){ return za/radius;}


// Dec 2 Prismatic& Prismatic::operator=(const Prismatic& pris)
// Dec 2 {
// Dec 2   std::cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<< std::endl;
// Dec 2   // super class
// Dec 2   this->m  =  pris.m;
// Dec 2   this->l  =  pris.l;
// Dec 2   this->r  =  pris.r;
// Dec 2   this->x0 = pris.x0;
// Dec 2   this->y0 = pris.y0;
// Dec 2   this->z0 = pris.z0;
// Dec 2   //	
// Dec 2   this->Iz = pris.Iz;
// Dec 2   this->radius = pris.radius;
// Dec 2   this->z  = pris.z ;
// Dec 2   this->zv = pris.zv;
// Dec 2   this->za = pris.za;
// Dec 2 
// Dec 2   std::cout << pris.m << "_" << (*this).m << std::endl;
// Dec 2   return *this;
// Dec 2 }

void Prismatic::setXYZ0(double _x0, double _y0, double _z0)
{
  x0 = _x0; y0 = _y0; z0 = _z0;
}

void Prismatic::setXYZ(double _x, double _y, double _z)
{
  //x = _x; y = _y; z = _z;
  z = _z;  
}

void Prismatic::setZsva(double s, double v, double a)
{
  z = s; zv = v; za = a; 
}


ScaraRobot::ScaraRobot(){}
ScaraRobot:: ScaraRobot(Revolute _r1, Revolute _r2, Prismatic _r3,
			Energy _e1, Energy _e2, Energy _e3)
  : r1(_r1), r2(_r2), r3(_r3), e1(_e1), e2(_e2), e3(_e3){}

ScaraRobot::ScaraRobot( const ScaraRobot& scara )
{
  *this = scara;
}

ScaraRobot::~ScaraRobot(){}

// int main(int argc, char *argv[])
// {
// 	Manipulator m1;
// 	Prismatic p1;
// 	Revolute r1;
// 	return 0;
// }

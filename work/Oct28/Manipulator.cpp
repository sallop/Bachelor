#include <iostream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include "Manipulator.hpp"

Manipulator::Manipulator()
{
  m  = 0.0; l  = 0.0; r  = 0.0;
  x0 = 0.0; y0 = 0.0; z0 = 0.0;
  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<<__LINE__<<std::endl;
  std::cout << "m=" << m << " l=" << l << " r=" << r << std::endl;
  std::cout << "x0=" << x0 << " y0=" << y0 << " z0=" << z0 << std::endl;
}

Manipulator::Manipulator(double _m, double _l,
			 double _x0, double _y0, double _z0)
{
  m  = _m; l   = _l; r   = _l / 2.0;
  x0 = _x0; y0 = _y0; z0 = _z0;
  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<<__LINE__<<std::endl;
  std::cout << "m=" << m << " l=" << l << " r=" << r << std::endl;
  std::cout << "x0=" << x0 << " y0=" << y0 << " z0=" << z0 << std::endl;
}

Manipulator::~Manipulator()
{
  std::cout << __FILE__<<"_"<<__FUNCTION__<<"_"<<__LINE__<<std::endl;
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
  std::cout <<__FILE__<<"_"<<__FUNCTION__<<"_"<<__LINE__<<std::endl;
}


VRManipulator::VRManipulator() : Manipulator()
{
  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__<< std::endl;
  x  = 0.0; xv  = 0.0; xa  = 0.0;
  y  = 0.0; yv  = 0.0; ya  = 0.0;
  th = 0.0; thv = 0.0; tha = 0.0;
  l  = 0.0; lv  = 0.0; la  = 0.0;
	
  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__ << std::endl;
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
	
  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__ << std::endl;
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

VRManipulator& VRManipulator::operator=(const VRManipulator& vrm)
{
  this->m  = vrm.m ;
  this->l  = vrm.l ;
  this->r  = vrm.r ;
  this->x0 = vrm.x0;
  this->y0 = vrm.y0;
  this->z0 = vrm.z0;
	
  this->x  = vrm.x ; this->xv  = vrm.xv ; this->xa  = vrm.xa ;
  this->y  = vrm.y ; this->yv  = vrm.yv ; this->ya  = vrm.ya ;
  this->th = vrm.th; this->thv = vrm.thv; this->tha = vrm.tha;
  this->l  = vrm.l ; this->lv  = vrm.lv ; this->la  = vrm.la ;
  //	this->z  = vrm.z ;
  return *this;
}

VRManipulator::~VRManipulator()
{
  std::cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<<std::endl;
}
				 

Revolute::Revolute() : Manipulator()
{
  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__ << std::endl;
}

Revolute::Revolute(double _m, double _l,
		   double _x0, double _y0, double _z0, double _rth)
  : Manipulator( _m,  _l,  _x0,  _y0,  _z0)
{
  Ig = (m*l*l) / 12.0;
  th = _rth; thv = 0.0;	tha = 0.0;
  x  =  0.0;  xv = 0.0;   xa = 0.0;
  y  =  0.0;  yv = 0.0;   ya = 0.0;
  std::cout << "Ig=" << Ig << std::endl;
  std::cout << "th=" << th << std::endl;
}

Revolute::Revolute(const Revolute &revo)
{
  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__ << std::endl;
  std::cout << "Copy constructor called." << std::endl;
  *this = revo;
}

Revolute& Revolute::operator=(const Revolute& revo)
{
  std::cout << __FILE__<<"_"<< __FUNCTION__ <<"_"<< __LINE__ <<std::endl;
  // super class
  this->m = revo.m;
  this->l = revo.l;
  this->r = revo.r;
  this->x0 = revo.x0;
  this->y0 = revo.y0;		
  this->z0 = revo.z0;	
  //
  this->Ig = revo.Ig;
	
  this->x  = revo.x ;
  this->xv = revo.xv;
  this->xa = revo.xa;
	
  this->y  = revo.y ;
  this->yv = revo.yv;
  this->ya = revo.ya;

  this->th  = revo.th ;
  this->thv = revo.thv;
  this->tha = revo.tha;
	
  return *this;
}

Revolute::~Revolute()
{
  std::cout << __FILE__ <<"_"<< __FUNCTION__ << "_" << __LINE__ << std::endl;
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
	
  std::cout << __FILE__     <<"_"
	    << __FUNCTION__ <<"_"
	    << __LINE__     << std::endl;
  std::cout << "x= " << x << std::endl;
  std::cout << "y= " << y << std::endl;
}

Prismatic::Prismatic() : Manipulator()
{
  std::cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<< std::endl;
}

Prismatic::Prismatic(double _m, double _l,
		     double _x0, double _y0, double _z0, double _radius)
  : Manipulator( _m,  _l,  _x0,  _y0,  _z0)
{
  radius = _radius;
  Iz = 0.0;
  z = z0 - l; zv = 0.0; za = 0.0;
		
  std::cout	<< "radius=" << radius << " Iz="<< Iz
		<< " z=" << z << " zv="<< zv << " za="<< za << std::endl;
}

Prismatic::Prismatic(const Prismatic &pris)
{
  std::cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<< std::endl;
  *this = pris;
}

Prismatic::~Prismatic()
{
  std::cout << __FILE__ <<"_"<< __FUNCTION__ <<"_" << __LINE__ << std::endl;
}

Prismatic& Prismatic::operator=(const Prismatic& pris)
{
  std::cout << __FILE__ <<"_"<<__FUNCTION__<<"_"<<__LINE__<< std::endl;
  // super class
  this->m  =  pris.m;
  this->l  =  pris.l;
  this->r  =  pris.r;
  this->x0 = pris.x0;
  this->y0 = pris.y0;
  this->z0 = pris.z0;
  //	
  this->Iz = pris.Iz;
  this->radius = pris.radius;
  this->z  = pris.z ;
  this->zv = pris.zv;
  this->za = pris.za;

  std::cout << pris.m << "_" << (*this).m << std::endl;
  return *this;
}

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


// int main(int argc, char *argv[])
// {
// 	Manipulator m1;
// 	Prismatic p1;
// 	Revolute r1;
// 	return 0;
// }

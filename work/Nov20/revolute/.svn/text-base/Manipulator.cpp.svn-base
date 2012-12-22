#include <iostream>
#include <cstdio>
#include <cmath>
#include "Manipulator.hpp"

using namespace std;

//cara::Scara(Manipulator *_r1,
//	     Manipulator *_r2,
//	     Manipulator *_r3)
// : r1(_r1), r2(_r2), r3(_r3)
//
//cout << __FILE__ <<" "<< __FUNCTION__ <<" "<< __LINE__ << endl;
//cout <<"r1="<<_r1 <<" "<<"r2="<<_r2 <<" "<<"r3="<<_r3 <<" "<< endl;  
//

ScaraTypeManipulator::ScaraTypeManipulator(Revolute  *_r1,
					   Revolute  *_r2,
					   Prismatic *_r3)
  : r1(_r1), r2(_r2), r3(_r3)
{
  cout << __FILE__ <<" "<< __FUNCTION__ <<" "<< __LINE__ << endl;
  cout <<"r1="<<_r1 <<" "<<"r2="<<_r2 <<" "<<"r3="<<_r3 <<" "<< endl;
}

//void Manipulator::print_value(FILE *_fp, char *name)
//{
//  fprintf(_fp, "variable:%s\n", name);
//  fprintf(_fp, "l=%lf,r=%lf,m=%lf,Ig=%lf\n", l, r, m, Ig);
//  fprintf(_fp, "x0=%lf,y0=%lf,z0=%lf\n", x0, y0, z0);
//  fprintf(_fp, "x  =%lf,y  =%lf,z  = %lf\n", x,  y,  z);
//  fprintf(_fp, "x1 =%lf,y1 =%lf,z1 = %lf\n", x1, y1, z1);
//  fprintf(_fp, "xv1=%lf,yv1=%lf\n",xv1, yv1);
//  fprintf(_fp, "xa1=%lf,ya1=%lf\n", xa1, ya1);
//  fprintf(_fp, "zv=%lf,za=%lf", zv, za);
//  fprintf(_fp, "th=%lf,thv=%lf,tha=%lf\n", th, thv, tha);
//  fprintf(_fp, "lv=%lf,la=%lf\n", lv, la);
//}

Manipulator::Manipulator()
  : m(0.0), r(0.0), Ig(0.0)
{
  cout << __FILE__ <<" "<< __FUNCTION__ <<" "<< __LINE__ << endl;  
}

Manipulator::Manipulator(double _m, double _l):
 m(_m),l(_l),r(_l/2.0),Ig(_m*_l*_l/12.0)
{
 cout << __FILE__ <<" "<< __FUNCTION__ <<" "<< __LINE__ << endl;
 cout << m <<" "<< l <<" "<< Ig <<" "<< endl;
}

void Manipulator::print_value(FILE *_fp, char *name)
{
  fprintf(_fp, "variable:%s\n", name);
  fprintf(_fp, "l=%lf,r=%lf,m=%lf,Ig=%lf\n", l, r, m, Ig);
  fprintf(_fp, "x0=%lf,y0=%lf,z0=%lf\n", x0, y0, z0);
  fprintf(_fp, "x  =%lf,y  =%lf,z  = %lf\n", x,  y,  z);
  fprintf(_fp, "x1 =%lf,y1 =%lf,z1 = %lf\n", x1, y1, z1);
  fprintf(_fp, "xv1=%lf,yv1=%lf\n",xv1, yv1);
  fprintf(_fp, "xa1=%lf,ya1=%lf\n", xa1, ya1);
  fprintf(_fp, "zv=%lf,za=%lf", zv, za);
  fprintf(_fp, "th=%lf,thv=%lf,tha=%lf\n", th, thv, tha);
  fprintf(_fp, "lv=%lf,la=%lf\n", lv, la);
}

void Manipulator::print_value(ostream& ofs, char *name)
{
  ofs << "variable:" << name << endl;
  ofs << "l="<<l<<",r="<<r<<",m="<<m<<",Ig="<<Ig << endl;
  ofs << "x0="<<x0<<",y0="<<y0<<",z0="<<z0<<endl;
  ofs << "x="<<x<<",y="<<y<<",z="<<z<<endl;
  ofs << "x1="<<x1<<",y1="<<y1<<",z1="<<z1<<endl;
  ofs << "xv1="<<xv1<<",yv1="<<yv1<<endl;
  ofs << "xa1="<<xa1<<",ya1="<<ya1<<endl;
  ofs << "zv="<<zv<<",za="<<za<<endl;
  ofs << "th="<<th<<",thv="<<thv<<",tha="<<tha<<endl;
  ofs << "lv="<<lv<<",la="<<la<<endl;
}

void Manipulator::setP0(double _x0, double _y0, double _z0)
{
  cout << __FILE__ <<" "<< __FUNCTION__ <<" "<< __LINE__ << endl;
  x0 = _x0; y0 = _y0; z0 = _z0;  
}

void Manipulator::setP1(double _AbsAngle){
  cout << __FILE__ <<" "<< __FUNCTION__ <<" "<< __LINE__ << endl;
  x1 = x0 + l*cos(_AbsAngle);
  y1 = y0 + l*sin(_AbsAngle);
  x = x0 + r*cos(_AbsAngle);
  y = y0 + r*sin(_AbsAngle);
  cout <<"x1="<<x1<<",y1="<<y1<<endl;
  cout <<"x ="<<x <<",y ="<<y <<endl;
}

void Manipulator::setTH(){
  cout << __FILE__ <<" "<< __FUNCTION__ <<" "<< __LINE__ << endl;  
}

AbstructManipulator::AbstructManipulator(double _m)
  : m(_m)
{
  cout << __FILE__ <<" "<< __FUNCTION__ <<" "<<__LINE__<< endl;
}

Revolute::Revolute(double _m, double _l)
  : AbstructManipulator(_m), l(_l), r(_l/2.0), Ig(_m*_l*_l/12.0)
{
  cout << __FILE__ <<" "<< __FUNCTION__ <<" "<<__LINE__<< endl;  
}

Prismatic::Prismatic(double _m, double _l)
  : AbstructManipulator(_m), l(_l), r(_l/2.0)
{
  cout << __FILE__ <<" "<< __FUNCTION__ <<" "<<__LINE__<< endl;    
}

//Manipulator man1(1.0, 1.0);
//int main(int argc, int *argv)
//{
//  Manipulator man2(2.0, 2.0);
//  Manipulator *man3 = new Manipulator(1.0, 1.0);
//  Scara sc(&man1, &man2, man3);
//  return 0;
//}

#ifndef MANIPULATOR_H
#define MANIPULATOR_H
#include <cstdio>
#include <iostream>

using namespace std;

struct Manipulator_;
struct Manipulator;
struct AbstructManipulator;
struct Prismatic;
struct Revolute;

struct Scara{
  Scara(Manipulator *r1, Manipulator *r2, Manipulator *r3);
  Manipulator *r1, *r2, *r3;
};

struct ScaraTypeManipulator{
  ScaraTypeManipulator(Revolute *r1, Revolute *r2, Prismatic *r3);
  Revolute  *r1, *r2;
  Prismatic *r3;
};

struct Manipulator_{
  double m;
  double Ig;
  double r;
  double l , lv, la;
  double x0, x, x1, xv1, xa1;
  double y0, y, y1, yv1, ya1;
  double z0, z, z1, zv , za ;
  double th, thv, tha;
  
  void print_value(FILE *_fp, char *name);
};

struct Manipulator{
  Manipulator();
  Manipulator(double m, double l);
//  const double m;
//  const double Ig;
//  const double r;
  double m;
  double Ig;
  double r;  
  double l , lv, la;
  double x0, x, x1, xv1, xa1;
  double y0, y, y1, yv1, ya1;
  double z0, z, z1, zv , za ;
  double th, thv, tha;
  
  void print_value(FILE *_fp, char *name);
  void print_value(ostream& ofs, char *name);
  void setP0(double x, double y, double z);// set point x0, y0
  void setP1(double absAngle);			// set point x1, y1
  void setTH();			// set relative angle th 
};

struct AbstructManipulator{
  AbstructManipulator(double m);
  const double m;
  virtual void print_value(FILE *_fp, char *name) = 0;
};

struct Revolute : AbstructManipulator{
  Revolute(double m, double l);
  const double Ig;
  const double r, l;

  double lv, la;
  double x0, x1, xv1, xa1;
  double y0, y1, yv1, ya1;
  double th, thv, tha;
  double z , z1, zv , za ;
};

struct Prismatic : AbstructManipulator{
  Prismatic(double m, double l);
  //  const double Ig;
  double l, r;
  double lv, la;
  double x;
  double y;
  double z0, z, z1, zv , za ;
  double th, thv, tha;		// using motor
};


#endif	// MANIPULATOR_H

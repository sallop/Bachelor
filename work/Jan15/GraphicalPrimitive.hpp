#ifndef __INCLUDE_GRAPHICALPRIMITIVE__
#define __INCLUDE_GRAPHICALPRIMITIVE__
#include <cstdio>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>

class GraphicalObject
{
public:
  virtual void draw(std::ostream& os) = 0;
  virtual ~GraphicalObject(){}
};

class Circle2D : public GraphicalObject
{
public:
  Circle2D(double _x=0, double _y=0, double _r=0);
  void setPoint(double x, double y);
  void setR(double r);
  void draw(std::ostream& os);
  double x, y, r;
private:
  double dstX[37], dstY[37];
};

class Circle3D : public GraphicalObject
{
public:
  Circle3D(double _x=0, double _y=0, double _z=0, double _r=0);
  void setPoint(double x, double y, double z);
  void setR(double r);
  void draw(std::ostream& os);
  double x, y, z, r;
private:
  double dstX[37], dstY[37];
};

class Line : public GraphicalObject
{
public:
  Line(double _x0, double _y0, double _x1, double _y1);
//   Line(double _x0, double _y0, double _th);  
  void setP0(double _x0, double _y0);
  void setP1(double _x1, double _y1);  
  void setP01(double _x0, double _y0, double _x1, double _y1);
  //  void setPoint(double _x0, double _y0, double _th);
  void draw(std::ostream& os);
private:
  double x0, y0, x1, y1, th;
};

class Line3D : public GraphicalObject
{
public:
  Line3D(double _x0, double _y0, double _z0,
	 double _x1, double _y1, double _z1);

  void setP0(double _x0, double _y0, double _z0);
  void setP1(double _x1, double _y1, double _z1);  
  void draw(std::ostream& os);
private:
  double x0, y0, z0, x1, y1, z1;
};

class Box : GraphicalObject{
public:
  Box(double x, double y, double r);
  void draw(std::ostream& os);
private:
  double x, y, r;
  double xs[4], ys[4];
};
#endif	//__INCLUDE_GRAPHICALPRIMITIVE__

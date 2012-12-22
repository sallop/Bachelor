#include "GraphicalPrimitive.hpp"
#define DEBUG 1

#if DEBUG
#include <fstream>
#endif

using namespace std; 

Circle2D::Circle2D(double _x, double _y, double _r)
  : x(_x), y(_y), r(_r)
{
  double th = 0.0;
  double dth = 2.0*M_PI/static_cast<double>(sizeof(dstX)/sizeof(*dstX));
  cout << __FILE__ << __FUNCTION__ << __LINE__ << endl;
  cout << x << y << r << endl;
  
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    dstX[i] = x + r*cos(th);
    dstY[i] = y + r*sin(th);
    th = th + dth;
  }
#if DEBUG
  cout <<"factor="<< sizeof(dstX)/sizeof(*dstX) << endl;
  FILE *fp = fopen("GraphicalTest.dat","w");
  fprintf(fp, "%lf %lf\n", x, y);  
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    fprintf(fp, "%lf %lf\n", dstX[i], dstY[i]);
  }
  fclose(fp);
#endif
  //  assert(false);
}

void Circle2D::setPoint(double _x, double _y)
{
  double dx = _x - x, dy = _y - y;
  cout <<"sizeof(dstX)="<< sizeof(dstX)
       <<"sizeof(*dstX)"<< sizeof(*dstX)
       << endl;
#if DEBUG
  FILE *fp = fopen("GraphicalTest.dat","w");
  fprintf(fp, "%lf %lf\n", x, y);
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    fprintf(fp,"%lf %lf\n", dstX[i], dstY[i]);
  }
#endif
  
  x = _x; y = _y;
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    dstX[i] = dstX[i] + dx;
    dstY[i] = dstY[i] + dy;
  }

#if DEBUG
  fprintf(fp,"%lf %lf\n", x, y);
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    fprintf(fp,"%lf %lf\n", dstX[i], dstY[i]);
  }
  fclose(fp);
#endif
  //  Assert(false);
}

void Circle2D::setR(double _r)
{
  double th = 0.0;// dth = 2.0*M_PI/36.0;
  double dth = 2.0*M_PI/(sizeof(dstX)/sizeof(*dstX));  
#if DEBUG
  FILE *fp=fopen("GraphicalTest.dat","w");
  for(int i=0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    fprintf(fp, "%lf %lf\n", dstX[i], dstY[i]);
  }
#endif

  r = _r;
  for(int i=0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    dstX[i] = x + r*cos(th);
    dstY[i] = y + r*sin(th);
    th = th + dth;
  }
  
#if DEBUG
  for(int i=0; i < (sizeof(dstX)/sizeof(*dstX)); ++i) {
    fprintf(fp, "%lf %lf\n", dstX[i], dstY[i]);
  }
  fclose(fp);
#endif
}

void Circle2D::draw(std::ostream &os)
{
  cout << x << setw(16) << y << endl;
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    os << dstX[i] << setw(16) << dstY[i] << endl;
  }
  os << dstX[0] << setw(16) << dstY[0] << endl << endl;  
}

Circle3D::Circle3D(double _x, double _y, double _z, double _r)
  : x(_x), y(_y), z(_z), r(_r)
{
  double  th = 0.0;
  double dth = 2.0*M_PI/static_cast<double>(sizeof(dstX)/sizeof(*dstX));
  cout << __FILE__ <<"_"<< __FUNCTION__ <<"_"<< __LINE__ << endl;
  cout << x <<"_"<< y <<"_" << r << endl;
  
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    dstX[i] = x + r*cos(th);
    dstY[i] = y + r*sin(th);
    th = th + dth;
  }
#if DEBUG
  cout <<"factor="<< sizeof(dstX)/sizeof(*dstX) << endl;
  FILE *fp = fopen("GraphicalTest.dat","w");
  fprintf(fp, "%lf %lf\n", x, y);  
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    fprintf(fp, "%lf %lf\n", dstX[i], dstY[i]);
  }
  fclose(fp);
#endif
  //  assert(false);
}

void Circle3D::setPoint(double _x, double _y, double _z)
{
  double dx = _x - x, dy = _y - y;
  cout <<"sizeof(dstX)="<< sizeof(dstX)
       <<"sizeof(*dstX)"<< sizeof(*dstX)
       << endl;
#if DEBUG
  FILE *fp = fopen("GraphicalTest.dat","w");
  fprintf(fp, "%lf %lf\n", x, y);
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    fprintf(fp,"%lf %lf\n", dstX[i], dstY[i]);
  }
#endif
  
  x = _x; y = _y;
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    dstX[i] = dstX[i] + dx;
    dstY[i] = dstY[i] + dy;
  }

#if DEBUG
  fprintf(fp,"%lf %lf\n", x, y);
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    fprintf(fp,"%lf %lf\n", dstX[i], dstY[i]);
  }
  fclose(fp);
#endif
  //  Assert(false);
}

void Circle3D::setR(double _r)
{
  double th = 0.0;// dth = 2.0*M_PI/36.0;
  double dth = 2.0*M_PI/(sizeof(dstX)/sizeof(*dstX));  
#if DEBUG
  FILE *fp=fopen("GraphicalTest.dat","w");
  for(int i=0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    fprintf(fp, "%lf %lf\n", dstX[i], dstY[i]);
  }
#endif

  r = _r;
  for(int i=0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    dstX[i] = x + r*cos(th);
    dstY[i] = y + r*sin(th);
    th = th + dth;
  }
  
#if DEBUG
  for(int i=0; i < (sizeof(dstX)/sizeof(*dstX)); ++i) {
    fprintf(fp, "%lf %lf\n", dstX[i], dstY[i]);
  }
  fclose(fp);
#endif
}

void Circle3D::draw(std::ostream &os)
{
  cout << x << setw(16) << y << endl;
  for(int i = 0; i < sizeof(dstX)/sizeof(*dstX); ++i){
    os << dstX[i] << setw(16)
       << dstY[i] << setw(16)
       <<    z    << std::endl;
  }
  os << dstX[0] << setw(16)
     << dstY[0] << setw(16)
     <<    z    << std::endl
     << std::endl;
}


Line::Line(double _x0, double _y0, double _x1, double _y1)
  : x0(_x0), y0(_y0), x1(_x1), y1(_y1)
{
  cout << __FUNCTION__ <<"-"<< __LINE__ << endl;
  cout <<x0<<" "<<y0<<" "<<x1<<" "<<y1<<endl;
}

void Line::setP0(double _x0, double _y0)
{
  x0 = _x0; y0 = _y0;
}

void Line::setP1(double _x1, double _y1)
{
  x1 = _x1; y1 = _y1;
}

void Line::setP01(double _x0, double _y0, double _x1, double _y1)
{
  x0 = _x0; y0 = _y0;
  x1 = _x1; y1 = _y1;
  cout <<x0<<" "<<y0<<" "<<x1<<" "<<" "<<y1<<endl;
}

void Line::draw(ostream &os)
{
  os << x0 << setw(16) << y0 << endl;
  os << x1 << setw(16) << y1 << endl;
  os << endl;
}

Line3D::Line3D(double _x0, double _y0, double _z0,
	       double _x1, double _y1, double _z1)
{
  x0 = _x0; y0 = _y0; z0 = _z0;
  x1 = _x1; y1 = _y1; z1 = _z1;
}

void Line3D::setP0(double _x0, double _y0, double _z0)
{
  x0 = _x0; y0 = _y0; z0 = _z0;  
}

void Line3D::setP1(double _x1, double _y1, double _z1)
{
  x1 = _x1; y1 = _y1; z1 = _z1;  
}

void Line3D::draw(std::ostream& os)
{
  os << x0 <<" "<< y0 <<" "<< z0 << std::endl;
  os << x1 <<" "<< y1 <<" "<< z1 << std::endl;
  os << std::endl;
}

Box::Box(double _x, double _y, double _r)
  : x(_x), y(_y), r(_r)
{
  double th = M_PI/4.0, dth = 2.0*M_PI/4.0;
  
  for(int i = 0; i < sizeof(xs)/sizeof(*xs); ++i){
    th = th + dth;
    xs[i] = x + r*cos(th);
    ys[i] = y + r*sin(th);
  }
}

void Box::draw(std::ostream& os)
{
  for(int i = 0; i < sizeof(xs)/sizeof(*xs); i++){
    os << xs[i] <<" "<< ys[i] << endl;
  }
  os << xs[0] <<" "<< ys[0] << endl << endl;  
}

 
// int main()
// {
//   double x = 0.0, y = 0.0;

//   Circle2D c0(x, y, 1.0);
//   ofstream ofs("test.dat");
//   c0.draw(ofs);
//   c0.setPoint(5.0, 1.0);
//   c0.draw(ofs);
//   c0.setR(2.0);
//   c0.draw(ofs);  
  
  
//   return 0;
// }

#ifndef BI3_2_H_
#define BI3_2_H_

class Bspline2d
{
public:
  double px[3], py[3], pz[3];
  double X(double t);
  double Y(double t);
  double Z(double t);
  Bspline2d();  
  Bspline2d(double x0, double y0, double z0,
	    double x1, double y1, double z1,
	    double x2, double y2, double z2 );
};
#endif	// BI3_2_H_

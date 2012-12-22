#include <iostream>
#include <cstdio>
#include "Bi3_2.hpp"

Bspline2d::Bspline2d(){}
Bspline2d::Bspline2d(double x0, double y0, double z0,
		     double x1, double y1, double z1,
		     double x2, double y2, double z2 )
{
  px[0] = x0; py[0] = y0; pz[0] = z0;
  px[1] = x1; py[1] = y1; pz[1] = z1;
  px[2] = x2; py[2] = y2; pz[2] = z2;
}

double Bspline2d::X(double t)
{
  return (1.0-t)*(1.0-t)*px[0] + 2.0*t*(1.0-t)*px[1] + t*t*px[2];
}

double Bspline2d::Y(double t)
{
  return (1.0-t)*(1.0-t)*py[0] + 2.0*t*(1.0-t)*py[1] + t*t*py[2];  
}

double Bspline2d::Z(double t)
{
  return (1.0-t)*(1.0-t)*pz[0] + 2.0*t*(1.0-t)*pz[1] + t*t*pz[2];  
}

// int main()
// {
//  double px[3] = {100, 200, 400};
//  double py[3] = {100, 500, 200};
//  double x, y; 
//  double step;
//  double t;

//  Bspline2d bspline2d(px[0], py[0],
// 		      px[1], py[1],
// 		      px[2], py[2]);
 
//  FILE *fp = fopen("dir-dat/Bi3_2.dat","w");
//  step = 1.0*0.01;
//  for(t = 0.0; t <= 1.0; t += step){
//    //    x = spline2d(t, px);
//    //    y = spline2d(t, py);
//    x = bspline2d.X(t);
//    y = bspline2d.Y(t);
//    fprintf(fp, "%lf %lf\n", x, y );
//  }
//  fclose(fp);

//  fp = popen("gnuplot", "w");
//  fprintf(fp, "set grid\n");  
//  fprintf(fp, "set size 0.6, 0.6\n");  
//  fprintf(fp, "plot 'dir-dat/Bi3_2.dat' u 1:2 w l\n");
//  fprintf(fp, "pause 1.5\n");

//  pclose(fp);

 
//  return 0;
//}

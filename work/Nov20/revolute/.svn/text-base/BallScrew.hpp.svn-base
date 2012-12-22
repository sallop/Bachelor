#ifndef BALLSCREW_H
#define BALLSCREW_H

class BallScrew{
public:
  BallScrew(double _J, double _m, double _pitch);
  double J;// motor side of inerita
  double m;// mass of fuka side
  double frac;// fraction of gear teeth  = 2.0*M_PI/pitch;
  double pitch;// pitch of screw
  void setJ(double _J);
  void setM(double _m);
  void setPitch(double _pitch);
  double getTau(double tha);// motor side of torque
};

#endif // BALLSCREW_H

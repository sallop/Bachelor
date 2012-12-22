#include <iostream>
#include <iomanip>
#include <cmath>
#include "BallScrew.hpp"

using namespace std;
BallScrew::BallScrew(double _J, double _m, double _pitch)
  : J(_J), m(_m), pitch(_pitch)
{
  frac = 2.0*M_PI/pitch;
  cout << "J=" << J << setw(2) 
       << "m=" << m << setw(2)
       << "f=" << frac << endl;
}

double BallScrew::getTau(double tha)
{
  return (J + m/(frac*frac))*tha;
}
void BallScrew::setJ(double _J){ J = _J;}
void BallScrew::setM(double _m){ m = _m;}
void BallScrew::setPitch(double _pitch)
{
  pitch = _pitch;
  frac = 2.0*M_PI/pitch;
}



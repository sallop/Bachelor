#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

class BallScrew{
public:
  BallScrew(double _J, double _m, double _pitch);
  double J;			// motor side of inerita
  double m;			// mass of fuka side
  double frac;			// fraction of gear teeth  = 2.0*M_PI/pitch;
  double pitch;			// pitch of screw
  double getTau(double tha);		// motor side of torque
};

BallScrew::BallScrew(double _J, double _m, double _pitch)
  : J(_J), m(_m), pitch(_pitch)
{
  frac = 2.0*M_PI/pitch;
  std::cout << "J=" << J << setw(2) 
	    << "m=" << m << setw(2)
	    << "f=" << frac << endl;
}

double BallScrew::getTau(double tha)
{
  return (J + m/(frac*frac))*tha;
}

class ReductionGearTrain{
public:
  ReductionGearTrain(double J1, double J2, double r1, double r2);  
  double J1;			// motor side of inerita
  double J2;			// fuka
  double r1;			// motor side of
  double r2;			//
  double frac;
  double getTorque(double tha);	// get torque of motor side
};

ReductionGearTrain::ReductionGearTrain(double _J1, double _J2,
				       double _r1, double _r2)
  : J1(_J1), J2(_J2), r1(_r1), r2(_r2), frac(r2/r1)
{
  cout <<"J1="<< J1 << " J2="<< J2 << " r1="<< r1 << " r2=" << r2 << endl;
}

double ReductionGearTrain::getTorque(double dth)
{
  return (J1 + J2/(frac*frac));
}

//class Mathstream : public ostream {
//  std::ostream& operator<< (std::ostream &os, double val){
//    return os << setw(4) << val << endl;
//  }
//};

int main(int argc, char **argv)
{
  std::ofstream ofs("BallScrew.dat");
  double dt, t;
  double dth, th;
  ReductionGearTrain rgt(1.0, 1.0, 1.0, 1.0);
  BallScrew bs(1.0, 1.0, 1.0);
  
  t = 0.0; th = 0.0; dt = 0.001; dth = M_PI/180.0;
  for(int i=0; i < 37; i++){
    t += dt;
    th += dth;
    ofs << t <<" "<< bs.getTau(i*dth) << endl;
  }
  
  return 0;
}

#ifndef __INCLUDE_ENERGY_H__
#define __INCLUDE_ENERGY_H__
#include <iostream>
#include <fstream>

struct Energy
{
  Energy();
  void plot(std::ofstream& ofs, double t);
  //  Energy& operator=(const Energy &ene);

  double ev, ia, tau, energy;
};
#endif // __INCLUDE_ENERGY_H__

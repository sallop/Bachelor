#include "Energy.hpp"
#include <cassert>

Energy::Energy()
  : ev(0.0), ia(0.0), tau(0.0), energy(0.0){}

void Energy::plot(std::ofstream& ofs, double t)
{
	ofs << t      <<"   " // 1
		<< ev     <<"   " // 2
		<< ia     <<"   " // 3
		<< tau    <<"   " // 4
		<< energy <<"   " // 5
		<< ev*ia  <<"   " // 6
		<< std::endl;
}

#include "Energy.hpp"
#include <cassert>

Energy::Energy()
{
  ev = 0.0;
  ia = 0.0;
  tau = 0.0;
  energy = 0.0;
}

Energy& Energy::operator=(const Energy& ene)
{
//printf("before\n");
//printf("%s_%s_%d\n", __FILE__, __FUNCTION__, __LINE__ );
//printf("this->\n");
//printf("ev=%lf,ia=%lf,tau=%lf,energy=%lf\n"
//	 , this->ev    
//	 , this->ia    
//	 , this->tau   
//	 , this->energy
//	 );
//printf("ene.\n");
//printf("ev=%lf,ia=%lf,tau=%lf,energy=%lf\n"
//	 , ene.ev    
//	 , ene.ia    
//	 , ene.tau   
//	 , ene.energy
//	 );
  this->ev    = ene.ev    ;
  this->ia    =	ene.ia    ;
  this->tau   =	ene.tau   ;
  this->energy=	ene.energy;

//  printf("after\n");
//  printf("%s_%s_%d\n", __FILE__, __FUNCTION__, __LINE__ );
//  printf("this->\n");
//  printf("ev=%lf,ia=%lf,tau=%lf,energy=%lf\n"
//	 , this->ev    
//	 , this->ia    
//	 , this->tau   
//	 , this->energy
//	 );
//  printf("ene.\n");
//  printf("ev=%lf,ia=%lf,tau=%lf,energy=%lf\n"
//	 , ene.ev    
//	 , ene.ia    
//	 , ene.tau   
//	 , ene.energy
//	 );
//  assert(false);
  return *this;
}

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

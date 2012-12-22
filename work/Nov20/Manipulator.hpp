#ifndef MANIPULATOR_H
#define MANIPULATOR_H
#include <cstdio>
#include <iostream>

using namespace std;

struct Manipulator;
struct Prismatic;
struct Revolute;

struct Manipulator
{
  Manipulator();
  Manipulator(double m, double l, double x0, double y0, double z0);
  virtual ~Manipulator();
  double m;		// mass of Manipulator; マニピュレータの質量
  double l, r;
  // long dimension of Manipulator; マニピュレータの長さ、　マニピュレータの長さ÷２
  double x0, y0, z0; // point of origin; 原点
  
  void print_value(std::ostream& os);
  //virtual void plot(std::ofstream& ofs);	//後で、各クラス毎に追加する
  virtual void setXYZ0(double _x0, double _y0, double _z0);
  virtual void setXYZ(double _x, double _y, double _z);
};

struct VRManipulator : Manipulator
{
  VRManipulator();
  VRManipulator(double _m, double _l,
		double _x0, double _y0, double _z0,
		double _rth);
  VRManipulator(const VRManipulator &vrm);
  ~VRManipulator();
  virtual void setXYZ0(double _x0, double _y0, double _z0);
  virtual void setXYZ(double _x, double _y, double _z);
  
  VRManipulator& operator=(const VRManipulator& vrm);
  double x , xv , xa ;
  double y , yv , ya ;
  double th, thv, tha;
  double l , lv , la ;
};

struct Revolute : Manipulator
{
  Revolute();
  Revolute(double m, double l, double x0, double y0, double z0, double rth);
  Revolute(const Revolute &revo); // copy constructor 初期化に用いられる
  ~Revolute();
  
  virtual void setXYZ0(double _x0, double _y0, double _z0);
  virtual void setXYZ(double _x, double _y, double _z);

  void setTHsva(double s, double v, double a);
  void setXsva(double s, double v, double a);
  void setYsva(double s, double v, double a);  
  
  Revolute& operator=(const Revolute&); // 代入演算子  式の途中での代入
  double Ig; // moment of inerita; 慣性モーメント
  double x , xv , xa ; // coordinate of endeffector; 先端、先端速度、先端加速度
  double y , yv , ya ; // coordinate of endeffector;
  double th, thv, tha;  // relative angle; ベースから先端までの相対角度
  
  void setEndEffector(double absTh); // absolute angle; ベースジョイントからの絶対角
};
 

struct Prismatic : Manipulator{
  Prismatic();
  Prismatic(double m, double l,
	    double x0, double y0, double z0,
	    double radius);
  Prismatic(const Prismatic &pris);	// copy constructor 初期化に使われる
  Prismatic& operator =(const Prismatic&); //　代入演算子のオーバーロード
  ~Prismatic();
  
  virtual void setXYZ0(double _x0, double _y0, double _z0);
  virtual void setXYZ(double _x, double _y, double _z);
  void setZsva(double s, double v, double a);
  
  double Iz; // moment of inerita; 慣性モーメント　(= 0)
  double radius; //　radius or pulley ;プーリーの半径
  double z, zv, za;
};
#endif	// MANIPULATOR_H

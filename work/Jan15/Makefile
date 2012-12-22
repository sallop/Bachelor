CC = g++
CFLAGS = -Wall -g -lm -O0

SRC = main.cpp Manipulator.cpp SinCurv.cpp Energy.cpp GraphicalPrimitive.cpp \
IDP.cpp IDPUP.cpp IDPUP_Mf.cpp IDPDOWN.cpp IDPDIRECT.cpp IDPNagato.cpp IDPZCurve.cpp variables.cpp
OBJ = $(SRC:%.cpp=%.o)

#OBJ = main.o Manipulator.o SinCurv.o Energy.o GraphicalPrimitive.o #ReadTrajectory.o

pulley : $(OBJ)
	$(CC) $^ -o $@ $(CFLAGS)
	echo "make $@ file. depends $^ files"

.cpp.o :
	echo "$< modified\n"
	$(CC) $< -c $(CFLAGS)

main.o : Manipulator.o SinCurv.o Energy.o variables.o IDP.o \
IDPUP.o IDPUP_Mf.o IDPDOWN.o IDPZCurve.o IDPNagato.o IDPDIRECT.o

variables.o : variables.cpp variables.hpp

IDP.o: IDP.hpp GraphicalPrimitive.hpp Manipulator.hpp Energy.hpp variables.hpp SinCurv.hpp
IDPZCurve.o: IDP.hpp GraphicalPrimitive.hpp Manipulator.hpp Energy.hpp variables.hpp SinCurv.hpp
IDPNagato.o : IDP.hpp GraphicalPrimitive.hpp Manipulator.hpp Energy.hpp variables.hpp SinCurv.hpp
IDPUP.o : IDP.hpp GraphicalPrimitive.hpp Manipulator.hpp Energy.hpp variables.hpp SinCurv.hpp
IDPUP_Mf.o : IDP.hpp GraphicalPrimitive.hpp Manipulator.hpp Energy.hpp variables.hpp SinCurv.hpp
IDPDOWN.o : IDP.hpp GraphicalPrimitive.hpp Manipulator.hpp Energy.hpp variables.hpp SinCurv.hpp
IDPDIRECT.o : IDP.hpp GraphicalPrimitive.hpp Manipulator.hpp Energy.hpp variables.hpp SinCurv.hpp

Manipulator.o : Manipulator.hpp
SinCurv.o : SinCurv.hpp
Energy.o :  Energy.hpp
GraphicalPrimitive.o :  GraphicalPrimitive.hpp

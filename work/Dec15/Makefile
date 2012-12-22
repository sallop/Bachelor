CC = g++
CFLAGS = -Wall -g -lm -O0
SRC = main.cpp Manipulator.cpp SinCurv.cpp Energy.cpp GraphicalPrimitive.cpp\
IDP.cpp variables.cpp
OBJ = $(SRC:%.c=%.o)

#OBJ = main.o Manipulator.o SinCurv.o Energy.o GraphicalPrimitive.o #ReadTrajectory.o

pulley : $(OBJ)
	$(CC) $^ -o $@ $(CFLAGS)
	echo "make $@ file. depends $^ files"

.cpp.o :
	echo "$< modified\n"
	$(CC) $< -c $(CFLAGS)

main.o : Manipulator.o SinCurv.o Energy.o variables.o IDP.o
variables.o : variables.cpp variables.hpp
IDP.o : IDP.cpp IDP.hpp Energy.hpp variables.hpp
Manipulator.o : Manipulator.cpp Manipulator.hpp
SinCurv.o : SinCurv.cpp SinCurv.hpp
Energy.o : Energy.cpp Energy.hpp
GraphicalPrimitive.o : GraphicalPrimitive.cpp GraphicalPrimitive.hpp



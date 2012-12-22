CC = g++
CFLAGS = -Wall -g -lm -O0
OBJ = main.o Manipulator.o SinCurv.o Energy.o GraphicalPrimitive.o

pulley : $(OBJ)
	$(CC) $^ -o $@ $(CFLAGS)

.cpp.o :
	echo "$< modified\n"
	$(CC) $< -c $(CFLAGS)

main.o : Manipulator.o SinCurv.o Energy.o
Manipulator.o : Manipulator.hpp
SinCurv.o : SinCurv.hpp
Energy.o : Energy.hpp
GraphicalPrimitive.o:GraphicalPrimitive.hpp
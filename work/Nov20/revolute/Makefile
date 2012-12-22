# Macro
PROG = revo
OBJ = main.o BallScrew.o Manipulator.o SinCurv.o

CC = g++
CFLAGS = -g -Wall

#
.SUFFIXES:.cpp.o
#
$(PROG):$(OBJ)
	$(CC) $(CFLAGS) -o $(PROG) $^
# suffix rule
.cpp.o:
	echo "making $<"
	$(CC) $(CFLAGS) -c $<

# remove target
.PHONY: clean
clean:
	ecno "remove $<"
	$(RM) $(PROG) $(OBJ)

# depend of header
main.o: BallScrew.o SinCurv.o
BallScrew.o: BallScrew.hpp
SinCurv.o: SinCurv.hpp

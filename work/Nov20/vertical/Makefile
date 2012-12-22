# define Macro
CC = g++
OBJS = main.o config.o Manipulator.o

# Make rule blocks
a.out : $(OBJS)
	echo "$? are modified."
	$(CC)  $(OBJS) -o $@

.cpp.o:
	$(CC) -c -g $<

main.o: config.o Manipulator.o
config.o: config.hpp
Manipulator.o: Manipulator.hpp

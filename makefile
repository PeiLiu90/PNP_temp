objects = PNP.o
CXX=g++
FC=gfortran
CFLAGS=-O3 

PNP.o:

main.o:main.cpp
	$(CXX) $(CFLAGS) -c main.cpp
run: $(objects) main.o
	$(CXX) $(CFLAGS) -o run $(objects) main.o

.PHONY:clean
clean:
	-rm run *.o *~



FC=ifort
FCFLAGS=-xHost -O2 -ipo -traceback -I/share/CGNS-intel/build/include -I../../src
#FCFLAGS=-g -fpe0 -traceback -I/share/CGNS-intel/build/include -I../../src
LIBNAME = libcns2d.a
LFLAGS=-L/share/CGNS-intel/build/lib -L../../lib
LIBS=-lcns2d -lcgns
TST_FILES = $(wildcard *.f90)
SRC_FILES = $(wildcard ../src/*.f90)
OBJ_FILES = $(SRC_FILES:.f90=.o)
TARGETS = ns_validation rectangle_grid

all : $(TARGETS)

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

ns_validation : ns_validation.o ../../lib/libcns2d.a
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

rectangle_grid : rectangle_grid.o
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

clean :
	rm -rf *.o *.out *.mod *.tec $(TARGETS)

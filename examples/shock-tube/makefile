FC=ifort
FCFLAGS=-xHost -O2 -ipo -traceback -I/share/CGNS-intel/build/include -I../../src
LIBNAME = libcns2d.a
LFLAGS=-L/share/CGNS-intel/build/lib -L../../lib
LIBS=-lcns2d -lcgns
TST_FILES = $(wildcard *.f90)
SRC_FILES = $(wildcard ../src/*.f90)
OBJ_FILES = $(SRC_FILES:.f90=.o)
TARGETS = shocktube shocktube_rotated tube_grid tube_grid_rotated

all : $(TARGETS)

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

shocktube : shocktube.o
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

shocktube_rotated : shocktube_rotated.o
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

tube_grid : tube_grid.o
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

tube_grid_rotated : tube_grid_rotated.o
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

clean :
	rm -rf *.o *.out *.mod *.tec $(TARGETS)

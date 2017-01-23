FC=ifort
FCFLAGS=-xHost -fp-model strict -O2 -ipo -traceback -I/share/CGNS-intel/build/include -I../../src
#FCFLAGS=-g -check all -fpe0 -warn -traceback -debug extended -I/share/CGNS-intel/build/include -I../src/
#FCFLAGS=-g -fpe0 -traceback -I/share/CGNS-intel/build/include -I../src/
LIBNAME = libcns2d.a
LFLAGS=-L/share/CGNS-intel/build/lib -L../../lib
LIBS=-lcns2d -lcgns
TST_FILES = $(wildcard *.f90)
SRC_FILES = $(wildcard ../src/*.f90)
OBJ_FILES = $(SRC_FILES:.f90=.o)
TARGETS = isentropic_vortex

all : $(TARGETS)

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

isentropic_vortex : isentropic_vortex.o
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

clean :
	rm -rf *.o *.out *.mod *.tec $(TARGETS)

FC=ifort
FCFLAGS=-O2 -I/share/CGNS-intel/build/include -I../src/ -mkl
LIBNAME = libeuler2d.a
LFLAGS=-L/share/CGNS-intel/build/lib -L../lib
LIBS=-leuler2d -lcgns
TST_FILES = $(wildcard *.f90)
SRC_FILES = $(wildcard ../src/*.f90)
OBJ_FILES = $(SRC_FILES:.f90=.o)
TARGETS = test_mesh test_cgns test_reflection test_ghost

all : lib $(TARGETS)

lib :
	cd ../src; make

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

test_cgns : test_cgns.o ../lib/$(LIBNAME)
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

test_mesh : test_mesh.o ../lib/$(LIBNAME)
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

test_reflection : test_reflection.o ../lib/$(LIBNAME)
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

test_ghost : test_ghost.o ../lib/$(LIBNAME)
	$(FC) $(FCFLAGS) -o $@ $< $(LFLAGS) $(LIBS)

clean :
	rm -rf *.o *.out *.mod *.tec $(TARGETS)

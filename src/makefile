FC=ifort
FCFLAGS=-xHost -O2 -ipo -traceback -I/share/CGNS-intel/build/include
#FCFLAGS=-g -fpe0 -traceback -I/share/CGNS-intel/build/include -I../src/
LIBNAME = libeuler2d
#SRC_FILES = $(wildcard ../src/*.f90)
SRC_FILES = utils.f90 mesh_class.f90 limiters.f90 grad.f90 flux.f90 riemann.f90 euler.f90
OBJ_FILES = $(SRC_FILES:.f90=.o)
TARGETS = $(LIBNAME)

all : $(TARGETS)

$(LIBNAME) : $(OBJ_FILES)
	ar rcs $@.a $^
	mkdir -p ../lib
	mv $@.a ../lib

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

clean :
	rm -rf *.o *.out *.mod

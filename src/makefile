FC=gfortran
#FCFLAGS=-g -check all -fpe0 -warn -traceback -debug extended -I/share/CGNS-intel/build/include -I../src/
FCFLAGS=-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -I/share/CGNS/build/include
LIBNAME = libeuler2d
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

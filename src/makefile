FC=ifort
FCFLAGS=-O2 -I/share/CGNS-intel/build/include -mkl
LIBNAME = libeuler2d
#SRC_FILES = $(wildcard ../src/*.f90)
SRC_FILES = utils.f90 mesh.f90 limiters.f90
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

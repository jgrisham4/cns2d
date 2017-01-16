FC=ifort
FCFLAGS=-fast -fp-model strict -traceback -I/share/CGNS-intel/build/include
#FCFLAGS=-g -fpe0 -traceback -I/share/CGNS-intel/build/include -I../src/ -check noarg_temp_created
LIBNAME = libcns2d
#SRC_FILES = $(wildcard ../src/*.f90)
LFLAGS = -L/share/CGNS-intel/build/lib -L../lib
LIBS = -lcns2d -lcgns
SRC_FILES = gas_properties.f90 mms.f90 utils.f90 mesh_class.f90 linalg.f90 solver_class.f90 flux.f90 riemann.f90 bcs.f90 acceleration.f90 limiters.f90 grad.f90 euler.f90 navierstokes.f90 temporal.f90
OBJ_FILES = $(SRC_FILES:.f90=.o)
TARGETS = $(LIBNAME) cns2d precns2d

all : $(TARGETS)

$(LIBNAME) : $(OBJ_FILES)
	ar rcs $@.a $^
	mkdir -p ../lib
	mv $@.a ../lib

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

cns2d : cns2d.o $(LIBNAME)
	$(FC) -o $@ $< $(FCFLAGS) $(LFLAGS) $(LIBS)

precns2d : precns2d.o $(LIBNAME)
	$(FC) -o $@ $< $(FCFLAGS) $(LFLAGS) $(LIBS)

clean :
	rm -rf *.o *.out *.mod cns2d precns2d
	rm -rf ../lib/*

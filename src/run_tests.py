#!/usr/bin/env python

import os
from subprocess import call

lapack95_path ="/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/include/intel64/lp64"

#funit_args = "FCFLAGS='-I/share/CGNS-intel/build/include -I{}' LDFLAGS='-L/share/CGNS-intel/build/lib -mkl'".format(lapack95_path)
#funit_args = "FCFLAGS='-I/share/CGNS-intel/build/include' LDFLAGS='-L/share/CGNS-intel/build/lib -mkl'"
funit_args = "FCFLAGS='-I/share/CGNS/build/include' LDFLAGS='-L/share/CGNS/build/lib'"
#funit_args = "FCFLAGS='-I/share/CGNS-intel/build/include' LDFLAGS='/share/CGNS-intel/build/lib/libcgns.a'"
#print("\ncalling funit with {}\n".format(funit_args))
call("source ~/.bashrc; source ~/.funit; {} funit".format(funit_args), shell=True)

#decision = input("\nWould you like to clean up the tests?\n(y/n): ")
decision = "y"
if decision == "y":
    call("source ~/.bashrc; source ~/.funit; funit --clean", shell=True)

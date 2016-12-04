#!/usr/bin/env python

import os
from subprocess import call

funit_args = "FCFLAGS='-I/share/CGNS-intel/build/include' LDFLAGS='-L/share/CGNS-intel/build/lib'"
#funit_args = "FCFLAGS='-I/share/CGNS-intel/build/include' LDFLAGS='/share/CGNS-intel/build/lib/libcgns.a'"
#print("\ncalling funit with {}\n".format(funit_args))
call("source ~/.bashrc; {} funit".format(funit_args), shell=True)

decision = input("\nWould you like to clean up the tests?\n(y/n): ")
if decision == "y":
    call("source ~/.bashrc; funit --clean", shell=True)

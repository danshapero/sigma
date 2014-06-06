SiGMA
=====

Sparse Graph and Matrix Algebra Library: the "i" is there to have a nice acronym. An old version of this library was called fempack before I came up with the snappy new name.

This is a library for performing operations on graphs and matrices. It is implemented
in Fortran 2003, with bindings for C and the accompanying family of languages.



Dependencies:
=============
* gfortran 4.7.2 or above, or ifort 13 or above; or, any Fortran compiler which can
  compile code written using the Fortran 2003 standard
* a C compiler
* OpenMP
* cmake


Installation:
=============
SiGMA is built using the program cmake. From the root SiGMA directory, the easiest
installation procedure is:
  mkdir build
  cd build
  cmake ..
  make
  make test
If you are compiling with the Intel Fortran compiler, make sure that your LD_LIBRARY_PATH variable is set to include the shared libraries provided with your installation of ifort.


Getting started:
================
In the directory examples/, you will find a series of sample programs which
demonstrate the capabilities and usage for many of the features of SiGMA. In order,
these are contained in the sub-directories:
  examples/graphs/
  examples/matrices/
  examples/solvers/
The user wishing to examine the implementations of the various features of the
library is advised to consult the corresponding folders in the src/ directory;
comments have been provided in order to illustrate the software design patterns
employed.


Digging deeper:
===============
The example programs serve to illustrate the basic functionalities of all the objects
defined in SiGMA. For the more advanced user, the programs in the folder apps/ exhibit
multiple features of the library.

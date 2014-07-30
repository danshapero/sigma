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
* mkdir build
* cd build
* cmake ..
* make
* make test
If you are compiling with the Intel Fortran compiler, make sure that your LD_LIBRARY_PATH variable is set to include the shared libraries provided with your installation of ifort.


Getting started:
================
In the directory examples/, you will find a series of sample programs which
demonstrate the capabilities and usage for many of the features of SiGMA. In order,
these are contained in the sub-directories:
* examples/graphs/
* examples/matrices/
* examples/solvers/
The user wishing to examine the implementations of the various features of the
library is advised to consult the corresponding folders in the src/ directory;
comments have been provided in order to illustrate the software design patterns
employed.


Digging deeper:
===============
The example programs serve to illustrate the basic functionalities of all the
objects defined in SiGMA. For the more advanced user, the programs in the folder
apps/ exhibit multiple features of the library.


Caveat emptor:
==============
SiGMA is still in development. Many features have yet to be implemented, and many others which have been implemented still need more guarantees on their correctness. The following things merit a warning:
* Never edit a graph while you're using a graph edge iterator. In most libraries that implement the iterator pattern for some container data type, they make sure that there's some exception which gets thrown or that the iterator is made invalid if you do this. SiGMA does not; it will simply fail, possibly unbeknownst to you.
* Ellpack graphs are a format especially suited to use on single-instruction, multiple-data (SIMD) architectures like graphics processing units or vector processors. However, many operations with ellpack graphs, e.g. edge iterators, fail if there is a vertex with no neighbors. If you use an ellpack graph, make sure there are no isolated vertices!
* There are numerous places where basic bounds checking could occur but does not. SiGMA will not check to make sure that all vectors, graphs and matrices are of compatible sizes, it will just fail or give you garbage results. If you want to test things, try using -DCMAKE_BUILD_TYPE=Debug when invoking cmake. This will set enable the -fbounds-check debug option, which catches a lot of bounds checking errors that could easily go unnoticed.
* There is as yet no routine to add a batch of edges to a graph, or to modify a batch of entries of a sparse matrix. Making it possible to execute lots of small operations instead as a big batch is crucial to cutting down on function call overhead.
* I got rid of parallel features during a refactor and haven't put them back in yet.
* There is no special consideration for symmetric matrices. While that may be advantageous for faster matrix multiplication, it also means that finding or modifying entries of the sparse matrix requires O(n) time as opposed to O(d), where d is the max degree of the underlying graph.
* Routines to compute the sum or product of two sparse matrices could be made faster by breaking the encapsulation of the sparse matrix object. Depending on how much of a bottleneck this operation proves to be, it may be worthwhile to overload those routines.
* Multiplying a row-oriented matrix with a column-oriented matrix is an extremely slow operation which requires a deep copy of one of the matrices into the opposite orientation. Never do it, ever.

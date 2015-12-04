This directory contains sources for operating on matrices, which extend the `linear_operator` interface with methods to modify matrix entries.


### Contents

* `sparse_matrix_interfaces.f90`: defines the `sparse_matrix_interface` class, which concrete implementations (CSR, ellpack, etc.) extend
* `formats/`: sources for concrete implementations of `sparse_matrix_interface`
* `sparse_matrix_factory.f90`: functions for instantiating sparse matrices in different formats
* `sparse_matrix_algebra.f90`: procedures for explicitly computing the sum or product of sparse matrices
Any other combination (row-row, col-col, col-row) is fine.
* `sparse_matrix_composites.f90`: defines the `sparse_matrix` class, which implements a block-storage format, i.e. a "matrix of matrices" 


### Design

The `sparse_matrix_interface` class serves as the [template class](https://en.wikipedia.org/wiki/Template_method_pattern) for the concrete implementations found in `formats/`.
Sparse matrices define iterators in a fashion similar to those of graphs, and a generic (if slow) implementation of matrix-vector multiplication can be defined in terms of the iterating over all the entries of the matrix.
The concrete matrix formats provide their own versions of matvec if they can do so faster than the default implementation.

The `sparse_matrix` class is a [composite](https://en.wikipedia.org/wiki/Composite_pattern) of sparse matrices.
This is particularly useful for problems that exhibit a natural coarse-grained block structure, e.g. the coupling of two distinct physical fields, or the coupling of a field with itself across an interface in the domain.
Moreover, the block structure can be nested.


### Usage

For best performance, a matrix should be built in two stages.
First, the connectivity graph `g` of the matrix is built in the same underlying format as the matrix (e.g. CS graph for a CSR/CSC matrix, etc.); and then the matrix `A` is constructed by invoking `A%set_graph(g)`.
This guarantees that there is no expensive intermediate stage of copying the graph structure from one format to another.
Moreover, if you have multiple matrices that share the same graph structure but with different entries, you can then share the graph `g` among all of them.

A matrix's entries can be changed after it has been built, but the user should take care to only set matrix entries `(i, j)` which are present in the connectivity graph; failure to do so will result in an expensive copy operation.

Avoid computing the product of a row-oriented matrix `A` and a column-oriented matrix `B`.
All other combinations (row/row, col/col, col/row) are fine.
When computing the conjugate product `transpose(P) * A * P` or `R * A * transpose(R)`, implemented in the routines `PtAP` and `RARt`, the appropriate multiplication order is selected in order to avoid this bad case.

!==========================================================================!
!==========================================================================!
module default_sparse_matrices                                             !
!==========================================================================!
!==========================================================================!
!====     This module contains the definition of the default           ====!
!==== implementation of sparse matrices. It uses no information about  ====!
!==== the underlying graph type other than what is provided by the     ====!
!==== graph interface. While it is very general, it is also not as     ====!
!==== fast as other, more specific implementations.                    ====!
!==========================================================================!
!==========================================================================!


use types, only: dp
use graph_interface
use generic_sparse_matrix_kernels
use sparse_matrix_interface

implicit none




!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: default_sparse_matrix                      !
!--------------------------------------------------------------------------!
    class(graph), pointer :: g
    real(dp), allocatable :: val(:)
contains


end type default_sparse_matrix



end module default_sparse_matrices

!==========================================================================!
!==========================================================================!
module sparse_matrices                                                     !
!==========================================================================!
!==========================================================================!

use graphs
use sparse_matrix_interface
use default_sparse_matrix_kernels
use default_matrices
use cs_matrices
use ellpack_matrices

implicit none

interface sparse_matrix
    module procedure sparse_matrix_factory
end interface


contains



!--------------------------------------------------------------------------!
function sparse_matrix_factory(nrow, ncol, g, orientation) result(A)       !
!--------------------------------------------------------------------------!
    integer, intent(in) :: nrow, ncol
    class(graph), pointer, intent(in) :: g
    character(len=3), intent(in) :: orientation
    class(sparse_matrix), pointer :: A

    select type(g)
        class is(cs_graph)
            A => cs_matrix(nrow, ncol, g, orientation)
        class is(ellpack_graph)
            A => ellpack_matrix(nrow, ncol, g, orientation)
        class default
            A => default_matrix(nrow, ncol, g, orientation)
    end select

end function sparse_matrix_factory





end module sparse_matrices

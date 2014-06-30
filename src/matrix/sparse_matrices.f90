!==========================================================================!
!==========================================================================!
module sparse_matrices                                                     !
!==========================================================================!
!==========================================================================!

use graphs
use sparse_matrix_interface
use default_sparse_matrix_kernels
use default_matrices

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

    allocate(default_matrix::A)

    select type(A)
        type is(default_matrix)
            call A%init(nrow, ncol, g, orientation)
    end select

end function sparse_matrix_factory





end module sparse_matrices

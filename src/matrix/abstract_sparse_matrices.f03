module abstract_sparse_matrices

use types
use abstract_graphs

implicit none



!--------------------------------------------------------------------------!
type, abstract :: sparse_matrix                                            !
!--------------------------------------------------------------------------!
    integer :: nrow, ncol, nnz
    logical :: pos_def
    class(graph), pointer :: g
contains
    procedure(matvec_ifc), deferred     :: matvec
end type sparse_matrix


abstract interface
    subroutine matvec_ifc(A,x,y)
        import :: sparse_matrix, dp
        class(sparse_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
    end subroutine matvec_ifc
end interface




end module abstract_sparse_matrices

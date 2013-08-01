module abstract_sparse_matrices

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
        import :: sparse_matrix
        class(sparse_matrix), intent(in) :: A
        real(kind(1d0)), intent(in)  :: x(:)
        real(kind(1d0)), intent(out) :: y(:)
    end subroutine matvec_ifc
end interface




end module abstract_sparse_matrices

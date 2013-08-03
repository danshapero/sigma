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
    logical :: assembled
contains
    procedure(init_ifc), deferred       :: init
    procedure(assemble_ifc), deferred   :: assemble
    procedure(get_value_ifc), deferred  :: get_value
    procedure(set_value_ifc), deferred  :: set_value, add_value
!    procedure(matrix_add_ifc), deferred :: matrix_add
    procedure(matvec_ifc), deferred     :: matvec
end type sparse_matrix


abstract interface
    subroutine init_ifc(A,nrow,ncol)
        import :: sparse_matrix
        class(sparse_matrix), intent(inout) :: A
        integer, intent(in) :: nrow, ncol
    end subroutine init_ifc

    subroutine assemble_ifc(A,g)
        import :: sparse_matrix, graph
        class(sparse_matrix), intent(inout) :: A
        class(graph), pointer, intent(in) :: g
    end subroutine assemble

    function get_value_ifc(A,i,j)
        import :: sparse_matrix, dp
        class(sparse_matrix), intent(in) :: A
        integer, intent(in) :: i,j
        real(dp) :: get_value_ifc
    end function get_value_ifc

    subroutine set_value_ifc(A,i,j,val)
        import :: sparse_matrix, dp
        class(sparse_matrix), intent(inout) :: A
        integer, intent(in) :: i,j
        real(dp), intent(in) :: val
    subroutine

    subroutine matvec_ifc(A,x,y)
        import :: sparse_matrix, dp
        class(sparse_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
    end subroutine matvec_ifc
end interface




end module abstract_sparse_matrices

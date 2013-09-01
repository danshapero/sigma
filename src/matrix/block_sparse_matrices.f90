module block_sparse_matrices

use types
use graphs
use sparse_matrices

implicit none


!--------------------------------------------------------------------------!
type, extends(sparse_matrix), abstract :: block_sparse_matrix              !
!--------------------------------------------------------------------------!
    integer :: nr, nc
contains
    procedure(get_block_ifc), deferred        :: get_block
    procedure(set_block_ifc), deferred        :: set_block
    procedure(set_block_ifc), deferred        :: add_block
    procedure(b_matvec_ifc), deferred         :: matvec
    procedure(b_matvec_ifc), deferred         :: matvec_t
    procedure(block_matvec_ifc), deferred     :: block_matvec
    procedure(block_matvec_ifc), deferred     :: block_matvec_t
    procedure(l_block_matvec_ifc), deferred   :: l_block_matvec
    procedure(l_block_matvec_ifc), deferred   :: l_block_matvec_t
    procedure(r_block_matvec_ifc), deferred   :: r_block_matvec
    procedure(r_block_matvec_ifc), deferred   :: r_block_matvec_t
!    generic :: matmul => matvec, &
!                 & block_matvec, & 
!               & l_block_matvec, &
!               & r_block_matvec
!    generic :: matmul_t => matvec_t, &
!                   & block_matvec_t, &
!                 & l_block_matvec_t, &
!                 & r_block_matvec_t
    generic :: matmul => block_matvec, &
         & l_block_matvec, &
         & r_block_matvec
    generic :: matmul_t => block_matvec_t, &
         & l_block_matvec_t, &
         & r_block_matvec_t
end type block_sparse_matrix


abstract interface
    subroutine init_block_matrix_ifc(A,nrow,ncol)
        import :: block_sparse_matrix
        class(block_sparse_matrix), intent(inout) :: A
        integer, intent(in) :: nrow, ncol
    end subroutine init_block_matrix_ifc

    subroutine block_assemble_ifc(A,g)
        import :: block_sparse_matrix, graph
        class(block_sparse_matrix), intent(inout) :: A
        class(graph), pointer, intent(in) :: g
    end subroutine block_assemble_ifc

    subroutine block_mat_neighbors_ifc(A,i,nbrs)
        import :: block_sparse_matrix
        class(block_sparse_matrix), intent(in) :: A
        integer, intent(in) :: i
        integer, intent(out) :: nbrs(:)
    end subroutine block_mat_neighbors_ifc

    function block_get_value_ifc(A,i,j)
        import :: block_sparse_matrix, dp
        class(block_sparse_matrix), intent(in) :: A
        integer, intent(in) :: i,j
        real(dp) :: block_get_value_ifc
    end function block_get_value_ifc

    subroutine block_set_value_ifc(A,i,j,val)
        import :: block_sparse_matrix, dp
        class(block_sparse_matrix), intent(inout) :: A
        integer, intent(in) :: i,j
        real(dp), intent(in) :: val
    end subroutine block_set_value_ifc

    subroutine b_matvec_ifc(A,x,y)
        import :: block_sparse_matrix, dp
        class(block_sparse_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
    end subroutine b_matvec_ifc

    subroutine get_block_ifc(A,i,j,vals)
        import :: block_sparse_matrix, dp
        class(block_sparse_matrix), intent(in) :: A
        integer, intent(in) :: i,j
        real(dp), intent(out) :: vals(:,:)
    end subroutine get_block_ifc

    subroutine set_block_ifc(A,i,j,vals)
        import :: block_sparse_matrix, dp
        class(block_sparse_matrix), intent(inout) :: A
        integer, intent(in) :: i,j
        real(dp), intent(in) :: vals(:,:)
    end subroutine set_block_ifc

    subroutine block_matvec_ifc(A,x,y)
        import :: block_sparse_matrix, dp
        class(block_sparse_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:,:)
        real(dp), intent(out) :: y(:,:)
    end subroutine block_matvec_ifc

    subroutine l_block_matvec_ifc(A,x,y)
        import :: block_sparse_matrix, dp
        class(block_sparse_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:,:)
    end subroutine l_block_matvec_ifc

    subroutine r_block_matvec_ifc(A,x,y)
        import :: block_sparse_matrix, dp
        class(block_sparse_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:,:)
        real(dp), intent(out) :: y(:)
    end subroutine r_block_matvec_ifc
end interface




end module block_sparse_matrices

!==========================================================================!
!==========================================================================!
module sparse_matrices                                                     !
!==========================================================================!
!==========================================================================!
!==== This module contains the definition of the basic real sparse     ====!
!==== matrix class. For more optimized formats, consider using block   ====!
!==== sparse matrices or variable-block sparse matrices.               ====!
!==========================================================================!
!==========================================================================!


use types, only: dp
use graphs
use ll_graphs

implicit none



!--------------------------------------------------------------------------!
type :: sparse_matrix                                                      !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    class(graph), pointer :: g
    integer :: nrow, ncol, nnz, max_degree, order(2)
    character(len=3) :: orientation
    logical :: pos_def
    logical :: assembled
contains
    procedure :: init => sparse_mat_init
!    Do we need this method at all?
!    procedure :: neighbors
    procedure :: get_value => sparse_mat_get_value
    procedure :: set_value => sparse_mat_set_value
    procedure :: add_value => sparse_mat_add_value
    procedure :: add_matrix => sparse_mat_add_mats
    procedure :: zero => sparse_mat_zero
    procedure :: left_permute => sparse_mat_leftperm
    procedure :: right_permute => sparse_mat_rightperm
    procedure :: matvec => sparse_mat_matvec
    procedure :: matvec_add => sparse_mat_matvec_add
    procedure :: destroy => sparse_mat_destroy
    generic :: matmul => matvec
end type sparse_matrix



!--------------------------------------------------------------------------!
type :: sparse_matrix_pointer                                              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), pointer :: A
end type sparse_matrix_pointer






contains



!--------------------------------------------------------------------------!
subroutine sparse_mat_init(A,nrow,ncol,orientation,g)
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol
    character(len=3), intent(in) :: orientation
    class(graph), pointer, intent(in), optional :: g

    ! Check if the user has provided a connectivity graph
    if (present(g)) then
        ! If so, make the matrix point to it
        A%g => g

        ! Check whether the matrix is stored in row- or column-major
        ! ordering
        select case(orientation)
            case('row')
                A%nrow = g%n
                A%ncol = g%m
            case('col')
                A%ncol = g%n
                A%nrow = g%m
        end select
    else
        ! If the user has not provided a connectivity graph, initialize the
        ! graph to a default type
        allocate(ll_graph::A%g)

        ! Set the number of left- and right-nodes for the graph to be the
        ! number of rows or columns of the matrix according to whether the
        ! matrix is row- or column-oriented
        select case(orientation)
            case('row')
                call A%g%init(nrow,ncol)
            case('col')
                call A%g%init(ncol,nrow)
        end select
    endif

    A%nnz = A%g%ne
    allocate(A%val(A%g%capacity))
    A%val = 0.0_dp
    A%max_degree = A%g%max_degree

    ! Set the order attribute for the matrix, so that the indices are
    ! reversed if we use column-major ordering
    ! I borrowed this approach from Numpy
    select case(orientation)
        case('row')
            A%order = [1, 2]
        case('col')
            A%order = [2, 1]
    end select

end subroutine sparse_mat_init



!--------------------------------------------------------------------------!
function sparse_mat_get_value(A,i,j) result(val)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: val
    ! local variables
    integer :: k, ind(2)

    ! (i,j) => (j,i) if the matrix is in column format
    ind = [i, j]
    ind = ind(A%order)

    ! Set the value to return to 0
    val = 0_dp

    ! Find the index k corresponding to the edge twixt (i,j) in A%g
    k = A%g%find_edge(ind(1),ind(2))

    ! If that edge exists, return the corresponding entry in the array
    ! of values of A
    if (k/=-1) val = A%val(k)

end function sparse_mat_get_value



!--------------------------------------------------------------------------!
subroutine sparse_mat_set_value(A,i,j,val)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k, ind(2)

    ! (i,j) => (j,i) if the matrix is in column format
    ind = [i, j]
    ind = ind(A%order)

    ! Find the index k corresponding to the edge twixt (i,j) in A%g
    k = A%g%find_edge(ind(1),ind(2))

    ! If that edge exists, set the corresponding entry in the array
    ! of values of A
    if (k/=-1) A%val(k) = val

end subroutine sparse_mat_set_value



!--------------------------------------------------------------------------!
subroutine sparse_mat_add_value(A,i,j,val)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k, ind(2)

    ! (i,j) => (j,i) if the matrix is in column format
    ind = [i, j]
    ind = ind(A%order)

    ! Find the index k corresponding to the edge twixt (i,j) in A%g
    k = A%g%find_edge(ind(1),ind(2))

    ! If that edge exists, set the corresponding entry in the array
    ! of values of A
    if (k/=-1) A%val(k) = A%val(k)+val

end subroutine sparse_mat_add_value



!--------------------------------------------------------------------------!
subroutine sparse_mat_add_mats(A,B,C)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    type(sparse_matrix), intent(in) :: B, C
    ! local variables
    integer :: i,j,k
    type(graph_edge_cursor) :: cursor

    ! Beats me

end subroutine sparse_mat_add_mats



!--------------------------------------------------------------------------!
subroutine sparse_mat_zero(A)                                              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A

    A%val = 0_dp

end subroutine sparse_mat_zero



!--------------------------------------------------------------------------!
subroutine sparse_mat_leftperm(A,p)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,j,k,ind(2)

    ! Uh...

end subroutine sparse_mat_leftperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_rightperm(A,p)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,j,k,ind(2)

end subroutine sparse_mat_rightperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_matvec(A,x,y,trans)                                  !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    logical, intent(in), optional :: trans

    y = 0_dp
    call A%matvec_add(x,y,trans)

end subroutine sparse_mat_matvec



!--------------------------------------------------------------------------!
subroutine sparse_mat_matvec_add(A,x,y,trans)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    logical, intent(in), optional :: trans
    ! local variables
    integer :: i,j,k,n,order(2),edges(2,64),num_returned,num_blocks
    type(graph_edge_cursor) :: cursor

    ! Reverse the matrix orientation if we're computing A'*x instead of A*x
    order = A%order
    if (present(trans)) then
        if (trans) order = [A%order(2),A%order(1)]
    endif

    ! set the parameters for the cursor
    cursor = A%g%make_cursor(0)

    ! Set the number of blocks in which we will be performing the matrix
    ! multiplication
    num_blocks = (cursor%final-cursor%start+1)/64+1

    do n=1,num_blocks
        ! Get the next 64 edges of the graph
        edges = A%g%get_edges(cursor,64,num_returned)

        ! Orient the edges according to the whether the matrix is in row-
        ! or column-orientation
        edges = edges(order,:)

        ! Go through all the edges (i,j) that we just plucked from A%g
        do k=1,num_returned
            i = edges(order(1),k)
            j = edges(order(2),k)

            ! Add A(i,j)*x(j) to y(i)
            if (i/=0 .and. j/=0) then
                y(i) = y(i)+A%val(64*(n-1)+k)*x(j)
            endif
        enddo
    enddo

end subroutine sparse_mat_matvec_add



!--------------------------------------------------------------------------!
subroutine sparse_mat_destroy(A)                                           !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A

    deallocate(A%val)
    nullify(A%g)

end subroutine sparse_mat_destroy




end module sparse_matrices

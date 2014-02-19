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
    procedure(sparse_mat_perm_ifc), pointer, private :: left_perm_impl
    procedure(sparse_mat_perm_ifc), pointer, private :: right_perm_impl
    procedure(sparse_matvec_add_ifc), pointer, private :: matvec_add_impl
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: init => sparse_mat_init


    !-----------
    ! Accessors
    !-----------
    procedure :: get_value => sparse_mat_get_value


    !----------
    ! Mutators
    !----------
    procedure :: set_value => sparse_mat_set_value
    procedure :: add_value => sparse_mat_add_value
    procedure :: add_matrix => sparse_mat_add_mats
    procedure :: zero => sparse_mat_zero
    procedure :: left_permute => sparse_mat_leftperm
    procedure :: right_permute => sparse_mat_rightperm
    procedure :: compress => sparse_mat_compress


    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure :: matvec => sparse_mat_matvec
    procedure :: matvec_add => sparse_mat_matvec_add


    !-------------
    ! Destructors
    !-------------
    procedure :: destroy => sparse_mat_destroy

end type sparse_matrix



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine sparse_mat_perm_ifc(A,p)
        import :: sparse_matrix
        class(sparse_matrix), intent(inout) :: A
        integer, intent(in) :: p(:)
    end subroutine sparse_mat_perm_ifc

    subroutine sparse_matvec_add_ifc(A,x,y,trans)
        import :: sparse_matrix, dp
        class(sparse_matrix), intent(in) :: A
        real(dp), intent(in)    :: x(:)
        real(dp), intent(inout) :: y(:)
        logical, intent(in), optional :: trans
    end subroutine sparse_matvec_add_ifc
end interface



!--------------------------------------------------------------------------!
type :: sparse_matrix_pointer                                              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), pointer :: A
end type sparse_matrix_pointer




contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

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
            A%left_perm_impl => sparse_mat_graph_leftperm
            A%right_perm_impl => sparse_mat_graph_rightperm
        case('col')
            A%order = [2, 1]
            A%left_perm_impl => sparse_mat_graph_rightperm
            A%right_perm_impl => sparse_mat_graph_leftperm
    end select

    ! Make the matrix's implementation of matvec refer to the decompressed
    ! version of the method
    A%matvec_add_impl => sparse_matvec_add_decompressed

end subroutine sparse_mat_init




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

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




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

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
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_perm_impl(p)

end subroutine sparse_mat_leftperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_rightperm(A,p)                                       !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_perm_impl(p)

end subroutine sparse_mat_rightperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_graph_leftperm(A,p)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i, source, dest, num_nodes
    integer, allocatable :: edge_p(:,:)
    real(dp), allocatable :: val(:)

    ! Permute the graph and get an array edge_p describing the permutation
    ! of the edges
    call A%g%left_permute(p,edge_p)

    ! If the edges require permutation,
    if (size(edge_p,2)>0) then
        ! Make a temporary array for the matrix values
        allocate(val(A%g%capacity))
        val = 0.0_dp

        ! Make an array of the re-arranged matrix values
        do i=1,size(edge_p,2)
            source = edge_p(1,i)
            dest = edge_p(2,i)
            num_nodes = edge_p(3,i)

            val(dest:dest+num_nodes-1) = A%val(source:source+num_nodes-1)
        enddo

        ! Transfer the allocation of the temporary array to the array
        ! of A's matrix entries
        call move_alloc(from=val, to=A%val)
    endif

    ! Deallocate the array describing the edge permutation
    deallocate(edge_p)

end subroutine sparse_mat_graph_leftperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_graph_rightperm(A,p)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i, source, dest, num_nodes
    integer, allocatable :: edge_p(:,:)
    real(dp), allocatable :: val(:)

    call A%g%right_permute(p,edge_p)

    if (size(edge_p,2)>0) then
        allocate(val(A%g%capacity))
        val = 0.0_dp

        do i=1,size(edge_p,2)
            source = edge_p(1,i)
            dest = edge_p(2,i)
            num_nodes = edge_p(3,i)

            val(dest:dest+num_nodes-1) = A%val(source:source+num_nodes-1)
        enddo

        call move_alloc(from=val, to=A%val)
    endif

    deallocate(edge_p)


end subroutine sparse_mat_graph_rightperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_compress(A)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    ! local variables
    integer :: i,source,dest,num
    integer, allocatable :: edge_p(:,:)
    real(dp), allocatable :: val(:)

    call A%g%compress(edge_p)

    if (size(edge_p,1)/=0) then
        ! Make a new array for the compressed matrix entries
        allocate(val(A%g%capacity))

        ! Move the uncompressed matrix entries into the compressed array
        do i=1,size(edge_p,1)
            source = edge_p(1,i)
            dest = edge_p(2,i)
            num = edge_p(3,i)

            val(dest:dest+num-1) = A%val(source:source+dest-1)
        enddo

        ! Transfer the allocation from the new array to the array that
        ! A points to
        call move_alloc(from=val, to=A%val)
    endif

    ! Change the strategy for the matrix's matvec to the implementation for
    ! compressed matrices, for which one needn't check
    A%matvec_add_impl => sparse_matvec_add_compressed

end subroutine sparse_mat_compress




!==========================================================================!
!==== Matrix-vector multiplication                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_mat_matvec(A,x,y,trans)                                  !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    logical, intent(in), optional :: trans

    y = 0_dp
    call A%matvec_add_impl(x,y,trans)

end subroutine sparse_mat_matvec



!--------------------------------------------------------------------------!
subroutine sparse_mat_matvec_add(A,x,y,trans)                              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    logical, intent(in), optional :: trans

    call A%matvec_add_impl(x,y,trans)

end subroutine sparse_mat_matvec_add



!--------------------------------------------------------------------------!
subroutine sparse_matvec_add_compressed(A,x,y,trans)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
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
        !edges = edges(order,:)

        ! Go through all the edges (i,j) that we just plucked from A%g
        do k=1,num_returned
            i = edges(order(1),k)
            j = edges(order(2),k)

            ! Add A(i,j)*x(j) to y(i). Since the underlying graph is
            ! compressed, we know that i,j are not zero.
            y(i) = y(i)+A%val(64*(n-1)+k)*x(j)
        enddo
    enddo

end subroutine sparse_matvec_add_compressed



!--------------------------------------------------------------------------!
subroutine sparse_matvec_add_decompressed(A,x,y,trans)                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
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
        !edges = edges(order,:)

        ! Go through all the edges (i,j) that we just plucked from A%g
        do k=1,num_returned
            i = edges(order(1),k)
            j = edges(order(2),k)

            ! Check to make sure that i,j are not zero, then add 
            ! A(i,j)*x(j) to y(i); a graph which has not yet been 
            ! compressed might have null edges
            if (i/=0 .and. j/=0) then
                y(i) = y(i)+A%val(64*(n-1)+k)*x(j)
            endif
        enddo
    enddo

end subroutine sparse_matvec_add_decompressed




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_mat_destroy(A)                                           !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A

    deallocate(A%val)
    nullify(A%g)

end subroutine sparse_mat_destroy




end module sparse_matrices

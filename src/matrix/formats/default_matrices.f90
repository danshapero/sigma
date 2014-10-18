!==========================================================================!
!==========================================================================!
module default_matrices                                                    !
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
use graph_interfaces
use ll_graphs
use default_sparse_matrix_kernels
use sparse_matrix_interfaces

implicit none




!--------------------------------------------------------------------------!
type, abstract, extends(sparse_matrix_interface) :: default_matrix         !
!--------------------------------------------------------------------------!
    class(graph_interface), pointer :: g => null()
    real(dp), allocatable :: val(:)
    integer :: ord(2)
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: copy_graph => default_matrix_copy_graph
    procedure :: set_graph => default_matrix_set_graph


    !-----------
    ! Accessors
    !-----------
    procedure :: get_nnz => default_matrix_get_nnz
    procedure :: get_value => default_matrix_get_value
    procedure :: get_row_degree => default_matrix_get_row_degree
    procedure :: get_column_degree => default_matrix_get_column_degree
    procedure :: get_row => default_matrix_get_row
    procedure :: get_column => default_matrix_get_column


    !-----------------------
    ! Edge, value iterators
    !-----------------------
    procedure :: make_cursor => default_matrix_make_cursor
    procedure :: get_edges => default_matrix_get_edges
    procedure :: get_entries => default_matrix_get_entries


    !----------
    ! Mutators
    !----------
    procedure :: set_value => default_matrix_set_value
    procedure :: add_value => default_matrix_add_value
    procedure :: zero => default_matrix_zero
    procedure :: left_permute => default_matrix_left_permute
    procedure :: right_permute => default_matrix_right_permute


    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    ! These procedures are provided by the parent class sparse_matrix.


    !-------------
    ! Destructors
    !-------------
    procedure :: destroy => default_matrix_destroy


    !--------------------
    ! Auxiliary routines
    !--------------------
    procedure(default_matrix_set_ordering_ifc), deferred :: set_ordering


    !--------------------------------------
    ! Implementations of matrix operations
    !--------------------------------------
    procedure(get_degree_kernel), deferred, nopass :: row_deg_impl
    procedure(get_degree_kernel), deferred, nopass :: col_deg_impl

    procedure(get_slice_kernel), deferred, nopass :: get_row_impl
    procedure(get_slice_kernel), deferred, nopass :: get_col_impl

    procedure(permute_kernel), deferred, nopass :: left_permute_impl
    procedure(permute_kernel), deferred, nopass :: right_permute_impl
end type default_matrix



!--------------------------------------------------------------------------!
type, extends(default_matrix) :: default_row_matrix                        !
!--------------------------------------------------------------------------!
contains
    !-----------
    ! Accessors
    !-----------
    procedure, nopass :: row_deg_impl => get_degree_contiguous
    procedure, nopass :: col_deg_impl => get_degree_discontiguous
    procedure, nopass :: get_row_impl => get_slice_contiguous
    procedure, nopass :: get_col_impl => get_slice_discontiguous

    !----------
    ! Mutators
    !----------
    procedure, nopass :: left_permute_impl  => graph_leftperm
    procedure, nopass :: right_permute_impl => graph_rightperm

    !--------------------
    ! Auxiliary routines
    !--------------------
    procedure :: set_ordering => default_row_matrix_set_ordering
end type default_row_matrix



!--------------------------------------------------------------------------!
type, extends(default_matrix) :: default_column_matrix                     !
!--------------------------------------------------------------------------!
contains
    !-----------
    ! Accessors
    !-----------
    procedure, nopass :: row_deg_impl => get_degree_discontiguous
    procedure, nopass :: col_deg_impl => get_degree_contiguous
    procedure, nopass :: get_row_impl => get_slice_discontiguous
    procedure, nopass :: get_col_impl => get_slice_contiguous

    !----------
    ! Mutators
    !----------
    procedure, nopass :: left_permute_impl  => graph_rightperm
    procedure, nopass :: right_permute_impl => graph_leftperm

    !--------------------
    ! Auxiliary routines
    !--------------------
    procedure :: set_ordering => default_column_matrix_set_ordering
end type default_column_matrix




!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine default_matrix_set_ordering_ifc(A)
        import :: default_matrix
        class(default_matrix), intent(inout) :: A
    end subroutine default_matrix_set_ordering_ifc
end interface


private :: default_matrix



contains





!==========================================================================!
!========                                                          ========!
!====         General routines for default row / column matrices       ====!
!========                                                          ========!
!==========================================================================!


!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine default_matrix_copy_graph(A, g)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(default_matrix), intent(inout) :: A
    class(graph_interface), intent(in) :: g
    ! local variables
    integer :: nv(2)
    logical :: trans

    call A%set_ordering()

    nv = [g%n, g%m]
    nv = nv(A%ord)

    if (A%nrow /= nv(1) .or. A%ncol /= nv(2)) then
        print *, 'Attemped to set default matrix connectivity structure to'
        print *, 'graph with inconsistent dimensions.'
        print *, 'Dimensions of matrix:', A%nrow, A%ncol
        print *, 'Dimensions of graph: ', nv(1), nv(2)
        call exit(1)
    endif

    ! Default to a list-of-lists graph for the underlying matrix format.
    ! We could choose the same format as the input graph, but this will
    ! work for the time being.
    !TODO: make it allocate A%g as a mold of g
    if (.not. associated(A%g)) allocate(ll_graph :: A%g)

    trans = A%ord(1) == 2
    call A%g%copy(g, trans)

    allocate(A%val(g%get_num_edges()))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine default_matrix_copy_graph



!--------------------------------------------------------------------------!
subroutine default_matrix_set_graph(A, g)                                  !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(inout) :: A
    class(graph_interface), target, intent(in) :: g

    call A%set_ordering()

    A%g => g
    call A%g%add_reference()

    allocate(A%val(g%get_num_edges()))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine default_matrix_set_graph




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function default_matrix_get_nnz(A) result(nnz)                             !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(in) :: A
    integer :: nnz

    nnz = A%g%get_num_edges()

end function default_matrix_get_nnz



!--------------------------------------------------------------------------!
function default_matrix_get_value(A, i, j) result(z)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(default_matrix), intent(in) :: A
    integer, intent(in) :: i, j
    real(dp) :: z
    ! local variables
    integer :: k, ind(2)

    ! Reverse the indices (i,j) according to the ordering (row- or column-
    ! major) of the matrix
    ind(1) = i
    ind(2) = j
    ind = ind(A%ord)

    ! Set the return value to 0
    z = 0.0_dp

    ! Find the index k in A%g of the edge twixt (i, j)
    k = A%g%find_edge(ind(1), ind(2))

    ! If that edge exists, find the corresponding matrix entry & return it
    if (k /= -1) z = A%val(k)

end function default_matrix_get_value



!--------------------------------------------------------------------------!
function default_matrix_get_row_degree(A, k) result(d)                     !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(in) :: A
    integer, intent(in) :: k
    integer :: d

    d = A%row_deg_impl(A%g, k)

end function default_matrix_get_row_degree



!--------------------------------------------------------------------------!
function default_matrix_get_column_degree(A, k) result(d)                  !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(in) :: A
    integer, intent(in) :: k
    integer :: d

    d = A%col_deg_impl(A%g, k)

end function default_matrix_get_column_degree



!--------------------------------------------------------------------------!
subroutine default_matrix_get_row(A, nodes, slice, k)                      !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k

    call A%get_row_impl( A%g, A%val, nodes, slice, k)

end subroutine default_matrix_get_row



!--------------------------------------------------------------------------!
subroutine default_matrix_get_column(A, nodes, slice, k)                   !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k

    call A%get_col_impl( A%g, A%val, nodes, slice, k)

end subroutine default_matrix_get_column




!==========================================================================!
!==== Edge, value iterators                                            ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function default_matrix_make_cursor(A) result(cursor)                      !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(in) :: A
    type(graph_edge_cursor) :: cursor

    cursor = A%g%make_cursor()

end function default_matrix_make_cursor



!--------------------------------------------------------------------------!
subroutine default_matrix_get_edges(A, edges, cursor, &                    !
                                                & num_edges, num_returned) !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(in) :: A
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned

    ! Get a batch of edges from the connectivity graph of A
    call A%g%get_edges(edges, cursor, num_edges, num_returned)

    ! Reverse the edges if A is in column-major order
    !TODO: See if we can do this without making a temporary array and if
    ! that's any faster
    edges = edges(A%ord, :)

end subroutine default_matrix_get_edges



!--------------------------------------------------------------------------!
subroutine default_matrix_get_entries(A, edges, entries, cursor, &         !
                                                & num_edges, num_returned) !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(default_matrix), intent(in) :: A
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    real(dp), intent(out) :: entries(num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned
    ! local variables
    integer :: indx

    ! Store the current position of the cursor
    indx = cursor%current

    ! Get the next batch of edges from A without the values
    call A%get_edges(edges, cursor, num_edges, num_returned)

    ! Get the entries from A
    entries = 0.0_dp
    entries(1 : num_returned) = A%val(indx + 1 : indx + num_returned)

end subroutine default_matrix_get_entries




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine default_matrix_set_value(A, i, j, z)                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(default_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k, ind(2)

    ind(1) = i
    ind(2) = j
    ind = ind(A%ord)

    k = A%g%find_edge(ind(1), ind(2))

    if (k /= -1) then
        A%val(k) = z
    else
        call set_matrix_value_with_reallocation(A%g, A%val, &
                                                    & ind(1), ind(2), z)
    endif

end subroutine default_matrix_set_value



!--------------------------------------------------------------------------!
subroutine default_matrix_add_value(A, i, j, z)                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(default_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k, ind(2)

    ind(1) = i
    ind(2) = j
    ind = ind(A%ord)

    k = A%g%find_edge(ind(1), ind(2))

    if (k /= -1) then
        A%val(k) = A%val(k) + z
    else
        call set_matrix_value_with_reallocation(A%g, A%val, &
                                                        & ind(1), ind(2), z)
    endif

end subroutine default_matrix_add_value



!--------------------------------------------------------------------------!
subroutine default_matrix_zero(A)                                          !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(inout) :: A

    A%val = 0.0_dp

end subroutine default_matrix_zero



!--------------------------------------------------------------------------!
subroutine default_matrix_left_permute(A, p)                               !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_permute_impl(A%g, A%val, p)

end subroutine default_matrix_left_permute



!--------------------------------------------------------------------------!
subroutine default_matrix_right_permute(A, p)                              !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_permute_impl(A%g, A%val, p)

end subroutine default_matrix_right_permute




!==========================================================================!
!==== Matrix-vector multiplication                                     ====!
!==========================================================================!

! The parent sparse matrix class provides a default, if not very optimal,
! implementation of these methods. The interface to the procedure is below
! but commented out, for the user who wishes to copy this module and use it
! as a guide to making her/his own matrix implementation.

!--------------------------------------------------------------------------!
! subroutine default_matrix_matvec_add(A, x, y)                            !
!--------------------------------------------------------------------------!
!     class(cs_matrix), intent(in) :: A                                    !
!     real(dp), intent(in)    :: x(:)                                      !
!     real(dp), intent(inout) :: y(:)                                      !
!                                                                          !
!     << Your implementation goes here >>                                  !
!                                                                          !
! end subroutine default_matrix_matvec_add                                 !
!--------------------------------------------------------------------------!



!--------------------------------------------------------------------------!
! subroutine default_matrix_matvec_t_add(A, x, y)                          !
!--------------------------------------------------------------------------!
!     class(cs_matrix), intent(in) :: A                                    !
!     real(dp), intent(in)    :: x(:)                                      !
!     real(dp), intent(inout) :: y(:)                                      !
!                                                                          !
!     << Your implementation goes here >>                                  !
!                                                                          !
! end subroutine default_matrix_matvec_t_add                               !
!--------------------------------------------------------------------------!




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine default_matrix_destroy(A)                                       !
!--------------------------------------------------------------------------!
    class(default_matrix), intent(inout) :: A

    ! Deallocate the array of A's matrix entries
    deallocate(A%val)

    ! Decrement the reference counter for A%g
    call A%g%remove_reference()

    ! Check how many references A%g has left -- if it's 0, we need to 
    ! destroy and deallocate it to avoid a memory leak
    if (A%g%reference_count == 0) then
        call A%g%destroy()
        deallocate(A%g)

    ! Otherwise, nullify the reference to the graph -- someone else is
    ! responsible for destroying it
    else
        nullify(A%g)
    endif

    A%graph_set = .false.
    A%dimensions_set = .false.

    A%reference_count = 0

end subroutine default_matrix_destroy





!==========================================================================!
!========                                                          ========!
!====    Format-specific routines for default row / column matrices    ====!
!========                                                          ========!
!==========================================================================!


!==========================================================================!
!==== Auxiliary routines                                               ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine default_row_matrix_set_ordering(A)                              !
!--------------------------------------------------------------------------!
    class(default_row_matrix), intent(inout) :: A

    A%ord = [1, 2]

end subroutine default_row_matrix_set_ordering



!--------------------------------------------------------------------------!
subroutine default_column_matrix_set_ordering(A)                           !
!--------------------------------------------------------------------------!
    class(default_column_matrix), intent(inout) :: A

    A%ord = [2, 1]

end subroutine default_column_matrix_set_ordering



end module default_matrices


!==========================================================================!
!==========================================================================!
module cs_matrices                                                         !
!==========================================================================!
!==========================================================================!
!====     This module contains the definition of compressed sparse     ====!
!==== matrices. These matrices explicitly use the compressed sparse    ====!
!==== graph data type, but perform matrix operations faster than the   ====!
!==== default matrix implementation.                                   ====!
!==========================================================================!
!==========================================================================!


use types, only: dp
use graph_interface
use cs_graphs
use default_sparse_matrix_kernels
use sparse_matrix_interface

implicit none




!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: cs_matrix                                  !
!--------------------------------------------------------------------------!
    class(cs_graph), pointer :: g
    real(dp), allocatable :: val(:)


    !---------------
    ! Function pointers for sundry matrix operations
    procedure(get_slice_kernel), pointer, nopass, private :: get_row_impl
    procedure(get_slice_kernel), pointer, nopass, private :: get_column_impl

    procedure(permute_kernel), pointer, nopass, private :: left_permute_impl
    procedure(permute_kernel), pointer, nopass, private :: right_permute_impl

    procedure(cs_matvec_add_ifc), pointer, private :: matvec_add_impl
    procedure(cs_matvec_add_ifc), pointer, private :: matvec_t_add_impl
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: set_ordering => cs_matrix_set_ordering
    procedure :: copy_graph_structure => cs_matrix_copy_graph_structure
    procedure :: set_graph => cs_matrix_set_graph


    !-----------
    ! Accessors
    !-----------
    procedure :: get_value => cs_matrix_get_value
    procedure :: get_row => cs_matrix_get_row
    procedure :: get_column => cs_matrix_get_column


    !-----------------------
    ! Edge, value iterators
    !-----------------------
    procedure :: make_cursor => cs_matrix_make_cursor
    procedure :: get_edges => cs_matrix_get_edges
    procedure :: get_entries => cs_matrix_get_entries


    !----------
    ! Mutators
    !----------
    procedure :: set_value => cs_matrix_set_value
    procedure :: add_value => cs_matrix_add_value
    procedure :: zero => cs_matrix_zero
    procedure :: left_permute => cs_matrix_left_permute
    procedure :: right_permute => cs_matrix_right_permute


    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure :: matvec_add   => cs_matvec_add
    procedure :: matvec_t_add => cs_matvec_t_add


    !-------------
    ! Destructors
    !-------------
    procedure :: destroy => cs_matrix_destroy

end type cs_matrix



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine cs_matvec_add_ifc(A, x, y)
        import :: cs_matrix, dp
        class(cs_matrix), intent(in) :: A
        real(dp), intent(in)    :: x(:)
        real(dp), intent(inout) :: y(:)
    end subroutine cs_matvec_add_ifc
end interface



contains





!==========================================================================!
!==== Constructors and factory methods                                 ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_matrix_set_ordering(A, orientation)                          !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    character(len=3), intent(in) :: orientation

    ! Set the `ord` attribute for the matrix, so that the indices are
    ! reversed if we use column-major ordering and so the right kernels
    ! are selected for matvec, getting rows/columns, etc.
    select case(orientation)
        case('row')
            A%ord = [1, 2]

            A%get_row_impl    => get_slice_contiguous
            A%get_column_impl => get_slice_discontiguous

            A%left_permute_impl  => graph_leftperm
            A%right_permute_impl => graph_rightperm

            A%matvec_add_impl   => csr_matvec_add
            A%matvec_t_add_impl => csc_matvec_add
        case('col')
            A%ord = [2, 1]

            A%get_row_impl    => get_slice_discontiguous
            A%get_column_impl => get_slice_contiguous

            A%left_permute_impl  => graph_rightperm
            A%right_permute_impl => graph_leftperm

            A%matvec_add_impl   => csc_matvec_add
            A%matvec_t_add_impl => csr_matvec_add
    end select

    A%ordering_set = .true.

end subroutine cs_matrix_set_ordering



!--------------------------------------------------------------------------!
subroutine cs_matrix_copy_graph_structure(A, g, trans)                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    class(graph), intent(in) :: g
    logical, intent(in), optional :: trans
    ! local variables
    integer :: ord(2), nv(2)
    logical :: tr

    tr = .false.
    if (present(trans)) tr = trans
    ord = [1, 2]
    if (tr) ord = [2, 1]

    nv = [g%n, g%m]
    nv = nv(ord)

    if (A%nrow /= nv(1) .or. A%ncol /= nv(2)) then
        print *, 'Attempted to set CS matrix connectivity structure to'
        print *, 'graph with inconsistent dimensions.'
        print *, 'Dimensions of matrix:', A%nrow, A%ncol
        print *, 'Dimensions of graph: ', nv(1), nv(2)
        call exit(1)
    endif

    allocate(A%g)
    call A%g%copy(g, tr)

    A%nnz = g%ne
    allocate(A%val(A%nnz))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine cs_matrix_copy_graph_structure



!--------------------------------------------------------------------------!
subroutine cs_matrix_set_graph(A, g)                                       !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    class(graph), target, intent(in) :: g

    if (A%nrow /= g%n .or. A%ncol /= g%m) then
        print *, 'Attempted to set CS matrix connectivity structure to'
        print *, 'graph with inconsistent dimensions.'
        print *, 'Dimensions of matrix:', A%nrow, A%ncol
        print *, 'Dimensions of graph: ', g%n, g%m
        call exit(1)
    endif

    select type(g)
        class is(cs_graph)
            A%g => g
        class default
            print *, 'Attempted to set CS matrix connectivity structure to'
            print *, 'point to a graph which is not a CS graph.'
            call exit(1)
    end select

    A%nnz = g%ne
    allocate(A%val(A%nnz))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine cs_matrix_set_graph




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function cs_matrix_get_value(A, i, j) result(z)                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
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

    do k = A%g%ptr(ind(1)), A%g%ptr(ind(1) + 1) - 1
        if ( A%g%node(k) == ind(2) ) z = A%val(k)
    enddo

end function cs_matrix_get_value



!--------------------------------------------------------------------------!
subroutine cs_matrix_get_row(A, nodes, slice, k)                           !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k

    call A%get_row_impl( A%g, A%val, nodes, slice, k)

end subroutine cs_matrix_get_row



!--------------------------------------------------------------------------!
subroutine cs_matrix_get_column(A, nodes, slice, k)                        !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k

    call A%get_column_impl( A%g, A%val, nodes, slice, k)

end subroutine cs_matrix_get_column




!==========================================================================!
!==== Edge, value iterators                                            ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function cs_matrix_make_cursor(A) result(cursor)                           !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    type(graph_edge_cursor) :: cursor

    cursor = A%g%make_cursor()

end function cs_matrix_make_cursor



!--------------------------------------------------------------------------!
subroutine cs_matrix_get_edges(A, edges, cursor, num_edges, num_returned)  !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned

    ! Get a batch of edges from the connectivity graph of A
    call A%g%get_edges(edges, cursor, num_edges, num_returned)

    ! Reverse the edges if A is in column-major order
    !TODO: See if we can do this without making a temporary array and if
    ! that's any faster
    edges = edges(A%ord, :)

end subroutine cs_matrix_get_edges



!--------------------------------------------------------------------------!
subroutine cs_matrix_get_entries(A, edges, entries, cursor, &              !
                                                & num_edges, num_returned) !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    integer, intent(out) :: edges(2, num_edges)
    real(dp), intent(out) :: entries(num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    ! local variables
    integer :: indx

    ! Store the current position of the cursor
    indx = cursor%current

    ! Get the next batch of eges from A without the values
    call A%get_edges(edges, cursor, num_edges, num_returned)

    ! Get the entries from A
    entries = 0.0_dp
    entries(1 : num_returned) = A%val(indx + 1 : indx + num_returned)

end subroutine cs_matrix_get_entries




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_matrix_set_value(A, i, j, z)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k, ind(2)
    logical :: found

    ind(1) = i
    ind(2) = j
    ind = ind(A%ord)

    found = .false.

    do k = A%g%ptr(ind(1)), A%g%ptr(ind(1) + 1) - 1
        if (A%g%node(k) == ind(2)) then
            A%val(k) = z
            found = .true.
        endif
    enddo

    if (.not. found) then
        call set_matrix_value_with_reallocation(A%g, A%val, &
                                                    & ind(1), ind(2), z)
        A%nnz = A%nnz + 1
    endif

end subroutine cs_matrix_set_value



!--------------------------------------------------------------------------!
subroutine cs_matrix_add_value(A, i, j, z)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k, ind(2)
    logical :: found

    ind(1) = i
    ind(2) = j
    ind = ind(A%ord)

    found = .false.

    do k = A%g%ptr(ind(1)), A%g%ptr(ind(1) + 1) - 1
        if (A%g%node(k) == ind(2)) then
            A%val(k) = A%val(k) + z
            found = .true.
        endif
    enddo

    if (.not. found) then
        call set_matrix_value_with_reallocation(A%g, A%val, &
                                                    & ind(1), ind(2), z)
        A%nnz = A%nnz + 1
    endif

end subroutine cs_matrix_add_value



!--------------------------------------------------------------------------!
subroutine cs_matrix_zero(A)                                               !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A

    A%val = 0.0_dp

end subroutine cs_matrix_zero



!--------------------------------------------------------------------------!
subroutine cs_matrix_left_permute(A, p)                                    !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_permute_impl(A%g, A%val, p)

end subroutine cs_matrix_left_permute



!--------------------------------------------------------------------------!
subroutine cs_matrix_right_permute(A, p)                                   !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_permute_impl(A%g, A%val, p)

end subroutine cs_matrix_right_permute




!==========================================================================!
!==== Matrix-vector multiplication                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_matvec_add(A, x, y)                                          !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)

    call A%matvec_add_impl(x, y)

end subroutine cs_matvec_add



!--------------------------------------------------------------------------!
subroutine cs_matvec_t_add(A, x, y)                                        !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)

    call A%matvec_t_add_impl(x, y)

end subroutine cs_matvec_t_add



!--------------------------------------------------------------------------!
subroutine csr_matvec_add(A, x, y)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i, j, k
    real(dp) :: z

    associate( g => A%g )


    do i = 1, g%n
        z = 0.0_dp

        do k = g%ptr(i), g%ptr(i + 1) - 1
            j = g%node(k)
            z = z + A%val(k) * x(j)
        enddo

        y(i) = y(i) + z
    enddo


    end associate

end subroutine csr_matvec_add



!--------------------------------------------------------------------------!
subroutine csc_matvec_add(A, x, y)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i, j, k
    real(dp) :: z

    associate( g => A%g )


    do j = 1, g%n
        z = x(j)

        do k = g%ptr(j), g%ptr(j + 1) - 1
            i = g%node(k)
            y(i) = y(i) + A%val(k) * z
        enddo
    enddo


    end associate

end subroutine csc_matvec_add




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_matrix_destroy(A)                                            !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A

    A%nnz = 0

    ! Deallocate the array of A's matrix entries
    deallocate(A%val)

    ! Decrement the reference counter for A%g
    call A%g%remove_reference()

    ! Nullify A's pointer to its graph. Don't de-allocate it -- there might
    ! still be other references to it someplace else.
    nullify(A%g)

end subroutine cs_matrix_destroy





end module cs_matrices

!==========================================================================!
!==========================================================================!
module ellpack_matrices                                                    !
!==========================================================================!
!==========================================================================!
!====     This module contains the definition of ellpack matrices.     ====!
!==== These matrices explicitly use the ellpack graph data type, but   ====!
!==== perform matrix operations faster than the default format.        ====!
!====     Ellpack matrices are especially suitable for execution on    ====!
!==== SIMD architectures due to the uniformity of the loops they use.  ====!
!==========================================================================!
!==========================================================================!


use types, only: dp
use graph_interface
use ellpack_graphs
use sparse_matrix_interface

implicit none




!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: ellpack_matrix                             !
!--------------------------------------------------------------------------!
    class(ellpack_graph), pointer :: g
    real(dp), allocatable :: val(:,:)


    !--------------------
    ! Function pointers for sundry matrix operations
    procedure(ellpack_get_slice_ifc), pointer, private :: get_row_impl
    procedure(ellpack_get_slice_ifc), pointer, private :: get_col_impl

    procedure(ellpack_permute_ifc), pointer, private :: left_permute_impl
    procedure(ellpack_permute_ifc), pointer, private :: right_permute_impl

    procedure(ellpack_matvec_add_ifc), pointer, private :: matvec_add_impl
    procedure(ellpack_matvec_add_ifC), pointer, private :: matvec_t_add_impl
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: init => ellpack_matrix_init


    !-----------
    ! Accessors
    !-----------
    procedure :: get_value => ellpack_matrix_get_value
    procedure :: get_row => ellpack_matrix_get_row
    procedure :: get_column => ellpack_matrix_get_column


    !-----------------------
    ! Edge, value iterators
    !-----------------------
    procedure :: make_cursor => ellpack_matrix_make_cursor
    procedure :: get_edges => ellpack_matrix_get_edges
    procedure :: get_entries => ellpack_matrix_get_entries


    !----------
    ! Mutators
    !----------
    procedure :: set_value => ellpack_matrix_set_value
    procedure :: add_value => ellpack_matrix_add_value
    procedure :: zero => ellpack_matrix_zero
    procedure :: left_permute => ellpack_matrix_left_permute
    procedure :: right_permute => ellpack_matrix_right_permute


    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure :: matvec_add   => ellpack_matvec_add
    procedure :: matvec_t_add => ellpack_matvec_t_add


    !-------------
    ! Destructors
    !-------------
    procedure :: destroy => ellpack_matrix_destroy


    !--------------------
    ! Auxiliary routines
    !--------------------
    procedure, private :: set_unallocated_matrix_value

end type ellpack_matrix



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine ellpack_get_slice_ifc(A, nodes, slice, k)
        import :: ellpack_matrix, dp
        class(ellpack_matrix), intent(in) :: A
        integer, intent(out) :: nodes(:)
        real(dp), intent(out) :: slice(:)
        integer, intent(in) :: k
    end subroutine ellpack_get_slice_ifc

    subroutine ellpack_permute_ifc(A, p)
        import :: ellpack_matrix
        class(ellpack_matrix), intent(inout) :: A
        integer, intent(in) :: p(:)
    end subroutine ellpack_permute_ifc

    subroutine ellpack_matvec_add_ifc(A, x, y)
        import :: ellpack_matrix, dp
        class(ellpack_matrix), intent(in) :: A
        real(dp), intent(in)    :: x(:)
        real(dp), intent(inout) :: y(:)
    end subroutine ellpack_matvec_add_ifc
end interface



interface ellpack_matrix
    module procedure ellpack_matrix_factory
end interface ellpack_matrix


private :: ellpack_matrix_leftperm, ellpack_matrix_rightperm


contains




!==========================================================================!
!==== Constructor and factory methods                                  ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_matrix_factory(nrow, ncol, g, orientation) result(A)      !
!--------------------------------------------------------------------------!
    integer, intent(in) :: nrow, ncol
    class(ellpack_graph), pointer, intent(in) :: g
    character(len=3), intent(in) :: orientation
    class(sparse_matrix), pointer :: A

    allocate(ellpack_matrix :: A)
    select type(A)
        class is(ellpack_matrix)
            call A%init(nrow, ncol, g, orientation)
    end select

end function ellpack_matrix_factory



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_init(A, nrow, ncol, g, orientation)              !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol
    class(ellpack_graph), pointer, intent(in) :: g
    character(len=3), intent(in) :: orientation

    A%nrow = nrow
    A%ncol = ncol

    A%g => g

    call A%g%add_reference()

    A%nnz = A%g%ne
    allocate(A%val(g%max_d, g%n))

    A%val = 0.0_dp

    ! Set the `ord` attribute for the matrix, so that the indices are
    ! reverse if we use column-major ordering and so the right kernels
    ! are selected for matvec, getting rows/columns, etc.
    select case(orientation)
        case('row')
            A%ord = [1, 2]

            A%get_row_impl => ellpack_get_slice_contiguous
            A%get_col_impl => ellpack_get_slice_discontiguous

            A%left_permute_impl  => ellpack_matrix_leftperm
            A%right_permute_impl => ellpack_matrix_rightperm

            A%matvec_add_impl   => ellpackr_matvec_add
            A%matvec_t_add_impl => ellpackc_matvec_add
        case('col')
            A%ord = [2, 1]

            A%get_row_impl => ellpack_get_slice_discontiguous
            A%get_col_impl => ellpack_get_slice_contiguous

            A%left_permute_impl  => ellpack_matrix_rightperm
            A%right_permute_impl => ellpack_matrix_leftperm

            A%matvec_add_impl   => ellpackc_matvec_add
            A%matvec_t_add_impl => ellpackr_matvec_add
    end select

end subroutine ellpack_matrix_init




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_matrix_get_value(A, i, j) result(z)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
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

    do k = 1, A%g%max_d
        if (A%g%node(k, ind(1)) == ind(2)) z = z + A%val(k, ind(1))
    enddo

end function ellpack_matrix_get_value



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_get_row(A, nodes, slice, k)                      !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k

    call A%get_row_impl(nodes, slice, k)

end subroutine ellpack_matrix_get_row



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_get_column(A, nodes, slice, k)                   !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k

    call A%get_col_impl(nodes, slice, k)

end subroutine ellpack_matrix_get_column



!--------------------------------------------------------------------------!
subroutine ellpack_get_slice_contiguous(A, nodes, slice, k)                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k
    ! local variables
    integer :: d

    nodes = 0
    slice = 0.0_dp

    d = A%g%degrees(k)

    nodes(1 : d) = A%g%node(1 : d, k)
    slice(1 : d) = A%val(1 : d, k)

end subroutine ellpack_get_slice_contiguous



!--------------------------------------------------------------------------!
subroutine ellpack_get_slice_discontiguous(A, nodes, slice, k)             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k
    ! local variables
    integer :: i, j, l, d, next

    ! Set the index in `nodes` and `slice` of the next entry to 0
    next = 0

    ! Set the returned nodes and values to 0
    nodes = 0
    slice = 0.0_dp

    ! Check every vertex `i`
    do i = 1, A%g%n
        ! Find the degree of `i`
        d = A%g%degrees(i)

        ! Check all the neighbors `j` of `i`
        do l = 1, d
            j = A%g%node(l, i)

            ! If `i` happens to neighbor `k`, put it and its matrix entry
            ! into the lists we're returning
            if (j == k) then
                next = next + 1

                ! Make `j` the next node of the slice
                nodes(next) = i

                ! and put in the corresponding matrix entry
                slice(next) = A%val(l, i)
            endif
        enddo
    enddo

end subroutine ellpack_get_slice_discontiguous




!==========================================================================!
!==== Edge, value iterators                                            ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_matrix_make_cursor(A) result(cursor)                      !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    type(graph_edge_cursor) :: cursor

    cursor = A%g%make_cursor()

end function ellpack_matrix_make_cursor



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_get_edges(A, edges, cursor, &                    !
                                                & num_edges, num_returned) !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned

    ! Get a batch of edges from the connectivity graph of `A`
    call A%g%get_edges(edges, cursor, num_edges, num_returned)

    ! Reverse the edges if `A` is in column-major order
    !TODO: see if we can do this without making a temporary array
    edges = edges(A%ord, :)

end subroutine ellpack_matrix_get_edges



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_get_entries(A, edges, entries, cursor, &         !
                                                & num_edges, num_returned) !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    integer, intent(out) :: edges(2, num_edges)
    real(dp), intent(out) :: entries(num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    ! local variables
    integer :: i, i1, i2, num_added, num_from_this_row, indx

    associate(g => A%g)

    !--------
    ! Get a batch of values from `A` and store them in the array `entries`

    entries = 0.0_dp

    num_returned = min(num_edges, cursor%final - cursor%current)

    indx = cursor%indx
    i1 = cursor%current / g%max_d + 1
    i2 = (cursor%current + num_returned - 1) / g%max_d + 1

    num_added = 0

    do i = i1, i2
        num_from_this_row = min( g%max_d - indx, num_returned - num_added)

        entries(num_added + 1 : num_added + num_from_this_row) &
                        & = A%val( indx + 1 : indx + num_from_this_row, i)
        num_added = num_added + num_from_this_row

        indx = mod(indx + num_from_this_row, g%max_d)
    enddo


    !--------
    ! Get a batch of edges from the connectivity graph of `A`

    call g%get_edges(edges, cursor, num_edges, num_returned)
    edges = edges(A%ord, :)

    end associate

end subroutine ellpack_matrix_get_entries




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_matrix_set_value(A, i, j, z)                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k, d, ind(2)
    logical :: found

    ind(1) = i
    ind(2) = j
    ind = ind(A%ord)

    found = .false.

    d = A%g%degrees(ind(1))
    do k = 1, d
        if (A%g%node(k, ind(1)) == ind(2)) then
            A%val(k, ind(1)) = z
            found = .true.
        endif
    enddo

    if (.not. found) call A%set_unallocated_matrix_value(ind(1), ind(2), z)

end subroutine ellpack_matrix_set_value



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_add_value(A, i, j, z)                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k, d, ind(2)
    logical :: found

    ind(1) = i
    ind(2) = j
    ind = ind(A%ord)

    found = .false.

    d = A%g%degrees(ind(1))
    do k = 1, d
        if (A%g%node(k, ind(1)) == ind(2)) then
            A%val(k, ind(1)) = A%val(k, ind(1)) + z
            found = .true.
        endif
    enddo

    if (.not. found) call A%set_unallocated_matrix_value(ind(1), ind(2), z)

end subroutine ellpack_matrix_add_value



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_zero(A)                                          !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A

    A%val = 0.0_dp

end subroutine ellpack_matrix_zero



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_left_permute(A, p)                               !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_permute_impl(p)

end subroutine ellpack_matrix_left_permute



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_right_permute(A, p)                              !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_permute_impl(p)

end subroutine ellpack_matrix_right_permute



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_leftperm(A, p)                                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i
    real(dp) :: val(A%g%max_d, A%g%n)

    val = A%val

    do i = 1, A%g%n
        A%val(:, p(i)) = val(:, i)
    enddo

    call A%g%left_permute(p)

end subroutine ellpack_matrix_leftperm



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_rightperm(A, p)                                  !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine ellpack_matrix_rightperm




!==========================================================================!
!==== Matrix-vector multiplication                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_matvec_add(A, x, y)                                     !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)

    call A%matvec_add_impl(x, y)

end subroutine ellpack_matvec_add



!--------------------------------------------------------------------------!
subroutine ellpack_matvec_t_add(A, x, y)                                   !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)

    call A%matvec_t_add_impl(x, y)

end subroutine ellpack_matvec_t_add



!--------------------------------------------------------------------------!
subroutine ellpackr_matvec_add(A, x, y)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i, j, k
    real(dp) :: z

    associate( g => A%g )


    do i = 1, g%n
        z = 0.0_dp

        do k = 1, g%max_d
            j = g%node(k, i)
            z = z + A%val(k, i) * x(j)            
        enddo

        y(i) = y(i) + z
    enddo


    end associate

end subroutine ellpackr_matvec_add



!--------------------------------------------------------------------------!
subroutine ellpackc_matvec_add(A, x, y)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i, j, k
    real(dp) :: z

    associate( g => A%g)


    do j = 1, g%n
        z = x(j)

        do k = 1, g%max_d
            i = g%node(k, j)
            y(i) = y(i) + A%val(k, j) * z
        enddo
    enddo


    end associate

end subroutine ellpackc_matvec_add




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_matrix_destroy(A)                                       !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A

    A%nnz = 0

    ! Deallocate the array of A's matrix entries
    deallocate(A%val)

    ! Decrement the reference counter for A%g
    call A%g%remove_reference()

    ! Nullify A's pointer to its graph. Don't de-allocate it -- there might
    ! still be other references to it someplace else.
    nullify(A%g)

end subroutine ellpack_matrix_destroy




!==========================================================================!
!==== Auxiliary routines                                               ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine set_unallocated_matrix_value(A, i, j, z)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: d
    real(dp), allocatable :: val(:, :)

    d = A%g%degrees(i)

    ! Expand the storage space for the values of `A` if need be
    if (d == A%g%max_d) then
        allocate(val(A%g%max_d + 1, A%g%n))
        val = 0.0_dp
        val(1 : A%g%max_d, :) = A%val(1 : A%g%max_d, :)
        call move_alloc(from = val, to = A%val)
    endif

    ! Add the edge to `g` and set the right entry in `val`
    call A%g%add_edge(i, j)
    A%val(d + 1, i) = z

end subroutine set_unallocated_matrix_value





end module ellpack_matrices

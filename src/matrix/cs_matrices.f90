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

    procedure(cs_matvec_add_kernel), pointer, private :: matvec_add_impl
    procedure(cs_matvec_add_kernel), pointer, private :: matvec_t_add_impl

contains
    !--------------
    ! Constructors
    !--------------
    procedure :: init => cs_matrix_init


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
    procedure :: assemble => cs_matrix_assemble
    procedure :: disassemble => cs_matrix_disassemble


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
    subroutine cs_matvec_add_kernel(A, x, y)
        import :: cs_matrix, dp
        class(cs_matrix), intent(in) :: A
        real(dp), intent(in)    :: x(:)
        real(dp), intent(inout) :: y(:)
    end subroutine cs_matvec_add_kernel
end interface



interface cs_matrix
    module procedure cs_matrix_factory
end interface




contains





!==========================================================================!
!==== Constructors and factory methods                                 ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function cs_matrix_factory(nrow, ncol, g, orientation) result(A)           !
!--------------------------------------------------------------------------!
    integer, intent(in) :: nrow, ncol
    class(cs_graph), pointer, intent(in) :: g
    character(len=3), intent(in) :: orientation
    class(sparse_matrix), pointer :: A

    allocate(cs_matrix :: A)
    select type(A)
        class is(cs_matrix)
            call A%init(nrow, ncol, g, orientation)
    end select

end function cs_matrix_factory



!--------------------------------------------------------------------------!
subroutine cs_matrix_init(A, nrow, ncol, g, orientation)                   !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol
    class(cs_graph), pointer, intent(in) :: g
    character(len=3), intent(in) :: orientation

    A%g => g

    call A%g%add_reference()

    A%nnz = A%g%ne
    allocate( A%val(A%g%capacity) )
    A%val = 0.0_dp

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
        case('col')
            A%ord = [2, 1]

            A%get_row_impl    => get_slice_discontiguous
            A%get_column_impl => get_slice_contiguous

            A%left_permute_impl  => graph_rightperm
            A%right_permute_impl => graph_leftperm
    end select

    A%matvec_add_impl   => sparse_matrix_matvec_add
    A%matvec_t_add_impl => sparse_matrix_matvec_t_add

end subroutine cs_matrix_init




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
    entries = A%val( indx + 1 : indx + num_returned )

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

    ind(1) = i
    ind(2) = j
    ind = ind(A%ord)

    do k = A%g%ptr(ind(1)), A%g%ptr(ind(1) + 1) - 1
        if (A%g%node(k) == ind(2)) A%val(k) = z
    enddo

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

    ind(1) = i
    ind(2) = j
    ind = ind(A%ord)

    do k = A%g%ptr(ind(1)), A%g%ptr(ind(1) + 1) - 1
        if (A%g%node(k) == ind(2)) A%val(k) = A%val(k) + z
    enddo

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



!--------------------------------------------------------------------------!
subroutine cs_matrix_assemble(A)                                           !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A

    ! Defer to a routine in default_sparse_matrix_kernels
    call assemble_matrix(A%g, A%val)

    ! Re-assign the matvec function pointers to the optimized routines
    select case(A%ord(1))
        case(1)
            A%matvec_add_impl   => csr_matvec_add
            A%matvec_t_add_impl => csc_matvec_add
        case(2)
            A%matvec_add_impl   => csc_matvec_add
            A%matvec_t_add_impl => csr_matvec_add
    end select

    A%assembled = .true.

end subroutine cs_matrix_assemble



!--------------------------------------------------------------------------!
subroutine cs_matrix_disassemble(A)                                        !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A

    ! Decompress the graph
    call A%g%decompress()

    ! Re-assign the matvec function pointers to the slow but safe routines
    A%matvec_add_impl   => sparse_matrix_matvec_add
    A%matvec_t_add_impl => sparse_matrix_matvec_t_add

    A%assembled = .false.

end subroutine cs_matrix_disassemble




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

    ! Deallocate the array of A's matrix entries
    deallocate(A%val)

    ! Decrement the reference counter for A%g
    call A%g%remove_reference()

    ! Nullify A's pointer to its graph. Don't de-allocate it -- there might
    ! still be other references to it someplace else.
    nullify(A%g)

end subroutine cs_matrix_destroy





end module cs_matrices

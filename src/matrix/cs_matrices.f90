!==========================================================================!
!==========================================================================!
module cs_matrices                                                         !
!==========================================================================!
!==========================================================================!
!====     This module contains the definition of compressed sparse     ====!
!==== matrices. These matrices explicitly use the compressed sparse    ====!
!==== graph data type, but perform matrix operations faster than the   ====!
!==== default matrix implementation.                                   ====!
!====     Compressed sparse matrices can be in row- or column-major    ====!
!==== ordering. These are modelled by the separate classes             ====!
!====         csr_matrix, csc_matrix.                                  ====!
!==== Many operations on CSR and CSC matrices do not depend on the     ====!
!==== matrix ordering; consequently, each type inherits from a common  ====!
!==== base class cs_matrix, which is private to this module.           ====!
!==========================================================================!
!==========================================================================!


use types, only: dp
use graph_interfaces
use cs_graphs
use default_sparse_matrix_kernels
use sparse_matrix_interfaces

implicit none




!--------------------------------------------------------------------------!
type, abstract, extends(sparse_matrix_interface) :: cs_matrix              !
!--------------------------------------------------------------------------!
    class(cs_graph), pointer :: g
    real(dp), allocatable :: val(:)
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: set_graph => cs_matrix_set_graph


    !-----------
    ! Accessors
    !-----------
    procedure :: get_row_degree => cs_matrix_get_row_degree
    procedure :: get_column_degree => cs_matrix_get_column_degree
    procedure :: get_row => cs_matrix_get_row
    procedure :: get_column => cs_matrix_get_column


    !-----------------------
    ! Edge, value iterators
    !-----------------------
    procedure :: make_cursor => cs_matrix_make_cursor
    procedure :: get_entries => cs_matrix_get_entries


    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure :: matvec_add => cs_matvec_add
    procedure :: matvec_t_add => cs_matvec_t_add


    !----------
    ! Mutators
    !----------
    procedure :: zero => cs_matrix_zero
    procedure :: left_permute => cs_matrix_left_permute
    procedure :: right_permute => cs_matrix_right_permute


    !-------------
    ! Destructors
    !-------------
    procedure :: destroy => cs_matrix_destroy


    !--------------------------------------
    ! Implementations of matrix operations
    !--------------------------------------
    procedure(cs_matvec_add_ifc), deferred, nopass :: matvec_add_impl
    procedure(cs_matvec_add_ifc), deferred, nopass :: matvec_t_add_impl

    procedure(get_degree_kernel), deferred, nopass :: row_deg_impl
    procedure(get_degree_kernel), deferred, nopass :: col_deg_impl

    procedure(get_slice_kernel), deferred, nopass :: get_row_impl
    procedure(get_slice_kernel), deferred, nopass :: get_col_impl

    procedure(permute_kernel), deferred, nopass :: left_permute_impl
    procedure(permute_kernel), deferred, nopass :: right_permute_impl

end type cs_matrix



!--------------------------------------------------------------------------!
type, extends(cs_matrix) :: csr_matrix                                     !
!--------------------------------------------------------------------------!
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: copy_graph_structure => csr_matrix_copy_graph_structure

    !-----------
    ! Accessors
    !-----------
    procedure :: get_value => csr_matrix_get_value
    procedure :: row_deg_impl => get_degree_contiguous
    procedure :: col_deg_impl => get_degree_discontiguous
    procedure :: get_row_impl => get_slice_contiguous
    procedure :: get_col_impl => get_slice_discontiguous

    !-----------------------
    ! Edge, value iterators
    !-----------------------
    procedure :: get_edges => csr_matrix_get_edges

    !----------
    ! Mutators
    !----------
    procedure :: set_value => csr_matrix_set_value
    procedure :: add_value => csr_matrix_add_value
    procedure :: left_permute_impl  => graph_leftperm
    procedure :: right_permute_impl => graph_rightperm

    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure :: matvec_add_impl   => csr_matvec_add
    procedure :: matvec_t_add_impl => csc_matvec_add
end type csr_matrix



!--------------------------------------------------------------------------!
type, extends(cs_matrix) :: csc_matrix                                     !
!--------------------------------------------------------------------------!
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: copy_graph_structure => csc_matrix_copy_graph_structure

    !-----------
    ! Accessors
    !-----------
    procedure :: get_value => csc_matrix_get_value
    procedure :: row_deg_impl => get_degree_discontiguous
    procedure :: col_deg_impl => get_degree_contiguous
    procedure :: get_row_impl => get_slice_discontiguous
    procedure :: get_col_impl => get_slice_contiguous

    !-----------------------
    ! Edge, value iterators
    !-----------------------
    procedure :: get_edges => csc_matrix_get_edges

    !----------
    ! Mutators
    !----------
    procedure :: set_value => csc_matrix_set_value
    procedure :: add_value => csc_matrix_add_value
    procedure :: left_permute_impl  => graph_rightperm
    procedure :: right_permute_impl => graph_leftperm

    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure :: matvec_add_impl   => csc_matvec_add
    procedure :: matvec_t_add_impl => csr_matvec_add
end type csc_matrix



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine cs_matvec_add_ifc(g, val, x, y)
        import :: cs_graph, dp
        class(cs_graph), intent(in) :: g
        real(dp), intent(in)    :: val(:)
        real(dp), intent(in)    :: x(:)
        real(dp), intent(inout) :: y(:)
    end subroutine cs_matvec_add_ifc
end interface


private :: cs_matrix



contains





!==========================================================================!
!==========================================================================!
!==== General routines for CS matrices                                 ====!
!==========================================================================!
!==========================================================================!


!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_matrix_set_graph(A, g)                                       !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    class(graph_interface), target, intent(in) :: g

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
function cs_matrix_get_row_degree(A, k) result(d)                          !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(in) :: k
    integer :: d

    d = A%row_deg_impl(A%g, k)

end function cs_matrix_get_row_degree



!--------------------------------------------------------------------------!
function cs_matrix_get_column_degree(A, k) result(d)                       !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(in) :: k
    integer :: d

    d = A%col_deg_impl(A%g, k)

end function cs_matrix_get_column_degree



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

    call A%get_col_impl( A%g, A%val, nodes, slice, k)

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

    call A%matvec_add_impl(A%g, A%val, x, y)

end subroutine cs_matvec_add



!--------------------------------------------------------------------------!
subroutine cs_matvec_t_add(A, x, y)                                        !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)

    call A%matvec_t_add_impl(A%g, A%val, x, y)

end subroutine cs_matvec_t_add




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




!==========================================================================!
!==========================================================================!
!==== Matrix-vector multiplication kernels                             ====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine csr_matvec_add(g, val, x, y)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    real(dp), intent(in)    :: val(:)
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i, j, k
    real(dp) :: z

    do i = 1, g%n
        z = 0.0_dp

        do k = g%ptr(i), g%ptr(i + 1) - 1
            j = g%node(k)
            z = z + val(k) * x(j)
        enddo

        y(i) = y(i) + z
    enddo

end subroutine csr_matvec_add



!--------------------------------------------------------------------------!
subroutine csc_matvec_add(g, val, x, y)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    real(dp), intent(in)    :: val(:)
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i, j, k
    real(dp) :: z

    do j = 1, g%n
        z = x(j)

        do k = g%ptr(j), g%ptr(j + 1) - 1
            i = g%node(k)
            y(i) = y(i) + val(k) * z
        enddo
    enddo

end subroutine csc_matvec_add




!==========================================================================!
!==========================================================================!
!==== CSR and CSC format-specific routines                             ====!
!==========================================================================!
!==========================================================================!


!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine csr_matrix_copy_graph_structure(A, g)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    class(graph), intent(in) :: g
    ! local variables

    if (A%nrow /= g%n .or. A%ncol /= g%m) then
        print *, 'Attempted to set CSR matrix connectivity structure to'
        print *, 'graph with inconsistent dimensions.'
        print *, 'Dimensions of matrix:', A%nrow, A%ncol
        print *, 'Dimensions of graph:', g%n, g%m
        call exit(1)
    endif

    if (.not. associated(A%g)) allocate(A%g)
    call A%g%copy(g)

    A%nnz = g%ne
    allocate(A%val(A%nnz))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine csr_matrix_copy_graph_structure



!--------------------------------------------------------------------------!
subroutine csc_matrix_copy_graph_structure(A, g)                           !
!--------------------------------------------------------------------------!
    class(csc_matrix), intent(inout) :: A
    class(graph), intent(in) :: g

    if (A%nrow /= g%m .or. A%ncol /= g%m) then
        print *, 'Attempted to set CSC matrix connectivity structure to'
        print *, 'graph with inconsistent dimensions.'
        print *, 'Dimensions of matrix:', A%nrow, A%ncol
        print *, 'Dimensions of graph:', g%m, g%n
        call exit(1)
    endif

    if (.not. associated(A%g)) allocate(A%g)
    call A%g%copy(g, trans = .true.)

    A%nnz = g%ne
    allocate(A%val(A%nnz))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine csc_matrix_copy_graph_structure




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function csr_matrix_get_value(A, i, j) result(z)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    integer, intent(in) :: i, j
    real(dp) :: z
    ! local variables
    integer :: k

    z = 0.0_dp

    do k = A%g%ptr(i), A%g%ptr(i + 1) - 1
        if (A%g%node(k) == j) z = A%val(k)
    enddo

end function csr_matrix_get_value



!--------------------------------------------------------------------------!
function csc_matrix_get_value(A, i, j) result(z)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    integer, intent(in) :: i, j
    real(dp) :: z
    ! local variables
    integer :: k

    z = 0.0_dp

    do k = A%g%ptr(j), A%g%ptr(j + 1) - 1
        if (A%g%node(k) == i) z = A%val(k)
    enddo

end function csc_matrix_get_value




!==========================================================================!
!==== Edge, value iterators                                            ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine csr_matrix_get_edges(A, edges, cursor, num_edges, num_returned) !
!--------------------------------------------------------------------------!
    class(csr_matrix), intent(in) :: A
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned

    call A%g%get_edges(edges, cursor, num_edges, num_returned)

end subroutine csr_matrix_get_edges



!--------------------------------------------------------------------------!
subroutine csc_matrix_get_edges(A, edges, cursor, num_edges, num_returned) !
!--------------------------------------------------------------------------!
    class(csc_matrix), intent(in) :: A
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned

    call A%g%get_edges(edges, cursor, num_edges, num_returned)
    edges = edges([2, 1], :)

end subroutine csc_matrix_get_edges




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine csr_matrix_set_value(A, i, j, z)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k
    logical :: found

    found = .false.

    do k = A%g%ptr(i), A%g%ptr(i + 1) - 1
        if (A%g%node(k) == j) then
            A%val(k) = z
            found = .true.
        endif
    enddo

    if (.not. found) then
        call set_matrix_value_with_reallocation(A%g, A%val, i, j, z)
        A%nnz = A%nnz + 1
    endif

end subroutine csr_matrix_set_value



!--------------------------------------------------------------------------!
subroutine csr_matrix_add_value(A, i, j, z)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k
    logical :: found

    found = .false.

    do k = A%g%ptr(i), A%g%ptr(i + 1) - 1
        if (A%g%node(k) == j) then
            A%val(k) = A%val(k) + z
            found = .true.
        endif
    enddo

    if (.not. found) then
        call set_matrix_value_with_reallocation(A%g, A%val, j, i, z)
        A%nnz = A%nnz + 1
    endif

end subroutine csr_matrix_add_value



!--------------------------------------------------------------------------!
subroutine csc_matrix_set_value(A, i, j, z)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csc_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k

    do k = A%g%ptr(j), A%g%ptr(j + 1) - 1
        if (A%g%node(k) == i) A%val(k) = z
    enddo

end subroutine csc_matrix_set_value



!--------------------------------------------------------------------------!
subroutine csc_matrix_add_value(A, i, j, z)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csc_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k

    do k = A%g%ptr(j), A%g%ptr(j + 1) - 1
        if (A%g%node(k) == i) A%val(k) = A%val(k) + z
    enddo

end subroutine csc_matrix_add_value



end module cs_matrices


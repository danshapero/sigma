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
use sparse_matrix_interfaces
use default_sparse_matrix_kernels

implicit none




!--------------------------------------------------------------------------!
type, abstract, extends(sparse_matrix_interface) :: cs_matrix              !
!--------------------------------------------------------------------------!
    class(cs_graph), pointer :: g => null()
    real(dp), allocatable :: val(:)
    logical, private :: get_row_is_fast = .false., get_col_is_fast = .false.
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: set_graph => cs_matrix_set_graph
    procedure :: copy_graph => cs_matrix_copy_graph


    !-----------
    ! Accessors
    !-----------
    procedure :: get_nnz => cs_matrix_get_nnz
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


    !--------------
    ! Optimization
    !--------------
    procedure :: is_get_row_fast => cs_matrix_is_get_row_fast
    procedure :: is_get_column_fast => cs_matrix_is_get_column_fast


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
    procedure :: set_dimensions => csr_matrix_set_dimensions

    !-----------
    ! Accessors
    !-----------
    procedure :: get_value => csr_matrix_get_value
    procedure, nopass :: row_deg_impl => get_degree_contiguous
    procedure, nopass :: col_deg_impl => get_degree_discontiguous
    procedure, nopass :: get_row_impl => get_slice_contiguous
    procedure, nopass :: get_col_impl => get_slice_discontiguous

    !-----------------------
    ! Edge, value iterators
    !-----------------------
    procedure :: get_edges => csr_matrix_get_edges

    !----------
    ! Mutators
    !----------
    procedure :: set_value => csr_matrix_set_value
    procedure :: add_value => csr_matrix_add_value
    procedure :: set_multiple_values => csr_matrix_set_multiple_values
    procedure :: add_multiple_values => csr_matrix_add_multiple_values
    procedure, nopass :: left_permute_impl  => graph_leftperm
    procedure, nopass :: right_permute_impl => graph_rightperm

    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure, nopass :: matvec_add_impl   => csr_matvec_add
    procedure, nopass :: matvec_t_add_impl => csc_matvec_add

end type csr_matrix



!--------------------------------------------------------------------------!
type, extends(cs_matrix) :: csc_matrix                                     !
!--------------------------------------------------------------------------!
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: set_dimensions => csc_matrix_set_dimensions

    !-----------
    ! Accessors
    !-----------
    procedure :: get_value => csc_matrix_get_value
    procedure, nopass :: row_deg_impl => get_degree_discontiguous
    procedure, nopass :: col_deg_impl => get_degree_contiguous
    procedure, nopass :: get_row_impl => get_slice_discontiguous
    procedure, nopass :: get_col_impl => get_slice_contiguous

    !-----------------------
    ! Edge, value iterators
    !-----------------------
    procedure :: get_edges => csc_matrix_get_edges

    !----------
    ! Mutators
    !----------
    procedure :: set_value => csc_matrix_set_value
    procedure :: add_value => csc_matrix_add_value
    procedure :: set_multiple_values => csc_matrix_set_multiple_values
    procedure :: add_multiple_values => csc_matrix_add_multiple_values
    procedure, nopass :: left_permute_impl  => graph_rightperm
    procedure, nopass :: right_permute_impl => graph_leftperm

    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure, nopass :: matvec_add_impl   => csc_matvec_add
    procedure, nopass :: matvec_t_add_impl => csr_matvec_add

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
!========                                                          ========!
!====             General routines for CSR / CSC matrices              ====!
!========                                                          ========!
!==========================================================================!


!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_matrix_copy_graph(A, g)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    class(graph_interface), intent(in) :: g
    ! local variables
    logical :: tr

    call check_source_dimensions(A, g%n, g%m)

    ! Allocate the matrix's connectivity structure if need be
    if (.not. associated(A%g)) allocate(A%g)

    ! Copy the input graph to the matrix's connectivity structure. If the
    ! matrix is in column-major ordering, we're actually copying the
    ! transpose of the graph.
    call A%g%copy(g, trans = A%get_col_is_fast)

    allocate(A%val(g%get_num_edges()))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine cs_matrix_copy_graph



!--------------------------------------------------------------------------!
subroutine cs_matrix_set_graph(A, g)                                       !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    class(graph_interface), target, intent(in) :: g

    select type(g)
        class is(cs_graph)
            A%g => g
            call A%g%add_reference()
        class default
            print *, 'Attempted to set CS matrix connectivity structure to'
            print *, 'point to a graph which is not a CS graph.'
            call exit(1)
    end select

    allocate(A%val(g%get_num_edges()))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine cs_matrix_set_graph




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function cs_matrix_get_nnz(A) result(nnz)                                  !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer :: nnz

    nnz = A%g%ne

end function cs_matrix_get_nnz



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
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    real(dp), intent(out) :: entries(num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned
    ! local variables
    integer :: idx

    ! Store the current position of the cursor
    idx = cursor%current

    ! Get the next batch of eges from A without the values
    call A%get_edges(edges, cursor, num_edges, num_returned)

    ! Get the entries from A
    entries = 0.0_dp
    entries(1 : num_returned) = A%val(idx+1 : idx+num_returned)

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

end subroutine cs_matrix_destroy




!==========================================================================!
!==== Optimization                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function cs_matrix_is_get_row_fast(A) result(fast)                         !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    logical :: fast

    fast = A%get_row_is_fast

end function cs_matrix_is_get_row_fast



!--------------------------------------------------------------------------!
function cs_matrix_is_get_column_fast(A) result(fast)                      !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    logical :: fast

    fast = A%get_col_is_fast

end function cs_matrix_is_get_column_fast





!==========================================================================!
!========                                                          ========!
!====               Matrix-vector multiplication kernels               ====!
!========                                                          ========!
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
!========                                                          ========!
!====               CSR / CSC format-specific routines                 ====!
!========                                                          ========!
!==========================================================================!


!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine csr_matrix_set_dimensions(A, nrow, ncol)                        !
!--------------------------------------------------------------------------!
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol

    A%nrow = nrow
    A%ncol = ncol

    A%dimensions_set = .true.

    call A%add_reference()

    A%get_row_is_fast = .true.
    A%get_col_is_fast = .false.

end subroutine csr_matrix_set_dimensions



!--------------------------------------------------------------------------!
subroutine csc_matrix_set_dimensions(A, nrow, ncol)                        !
!--------------------------------------------------------------------------!
    class(csc_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol

    A%nrow = nrow
    A%ncol = ncol

    A%dimensions_set = .true.

    call A%add_reference()

    A%get_row_is_fast = .false.
    A%get_col_is_fast = .true.

end subroutine csc_matrix_set_dimensions




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
    class(csc_matrix), intent(in) :: A
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
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned

    call A%g%get_edges(edges, cursor, num_edges, num_returned)

end subroutine csr_matrix_get_edges



!--------------------------------------------------------------------------!
subroutine csc_matrix_get_edges(A, edges, cursor, num_edges, num_returned) !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csc_matrix), intent(in) :: A
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned
    ! local variables
    integer :: j, k, num, bt

    associate(current => cursor%current, g => A%g)

    edges = 0
    num_returned = min(num_edges, cursor%last - current)

    edges(1, 1 : num_returned) = g%node(current+1 : current+num_returned)

    j = cursor%edge(1)

    k = 0
    do while(k < num_returned)
        num = min(g%ptr(j + 1) - 1 - current, num_returned - k)

        edges(2, k+1 : k+num) = j

        bt = (sign(1, num - (g%ptr(j+1) - 1 - current)) + 1) / 2
        j = j + bt

        current = current + num
        k = k + num
    enddo

    cursor%edge(1) = j

    end associate

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
        call set_matrix_value_with_reallocation(A%g, A%val, i, j, z)
    endif

end subroutine csr_matrix_add_value



!--------------------------------------------------------------------------!
subroutine csr_matrix_set_multiple_values(A, is, js, B)                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: is(:), js(:)
    real(dp), intent(in) :: B(:,:)
    ! local variables
    integer :: i, j, k, l, m
    real(dp) :: z
    logical :: found

    do k = 1, size(is)
        i = is(k)

        do l = 1, size(js)
            j = js(l)
            z = B(k, l)

            found = .false.

            do m = A%g%ptr(i), A%g%ptr(i + 1) - 1
                if (A%g%node(m) == j) then
                    A%val(m) = z
                    found = .true.
                endif
            enddo

            if (.not. found) then
                call set_matrix_value_with_reallocation(A%g, A%val, i, j, z)
            endif
        enddo
    enddo

end subroutine csr_matrix_set_multiple_values



!--------------------------------------------------------------------------!
subroutine csr_matrix_add_multiple_values(A, is, js, B)                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: is(:), js(:)
    real(dp), intent(in) :: B(:,:)
    ! local variables
    integer :: i, j, k, l, m
    real(dp) :: z
    logical :: found

    do k = 1, size(is)
        i = is(k)

        do l = 1, size(js)
            j = js(l)
            z = B(k, l)

            found = .false.

            do m = A%g%ptr(i), A%g%ptr(i + 1) - 1
                if (A%g%node(m) == j) then
                    A%val(m) = A%val(m) + z
                    found = .true.
                endif
            enddo

            if (.not. found) then
                call set_matrix_value_with_reallocation(A%g, A%val, i, j, z)
            endif
        enddo
    enddo

end subroutine csr_matrix_add_multiple_values



!--------------------------------------------------------------------------!
subroutine csc_matrix_set_value(A, i, j, z)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csc_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k
    logical :: found

    found = .false.

    do k = A%g%ptr(j), A%g%ptr(j + 1) - 1
        if (A%g%node(k) == i) then
            A%val(k) = z
            found = .true.
        endif
    enddo

    if (.not. found) then
        call set_matrix_value_with_reallocation(A%g, A%val, j, i, z)
    endif

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
    logical :: found

    found = .false.

    do k = A%g%ptr(j), A%g%ptr(j + 1) - 1
        if (A%g%node(k) == i) then
            A%val(k) = A%val(k) + z
            found = .true.
        endif
    enddo

    if (.not. found) then
        call set_matrix_value_with_reallocation(A%g, A%val, j, i, z)
    endif

end subroutine csc_matrix_add_value



!--------------------------------------------------------------------------!
subroutine csc_matrix_set_multiple_values(A, is, js, B)                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csc_matrix), intent(inout) :: A
    integer, intent(in) :: is(:), js(:)
    real(dp), intent(in) :: B(:,:)
    ! local variables
    integer :: i, j, k, l, m
    real(dp) :: z
    logical :: found

    do l = 1, size(js)
        j = js(l)

        do k = 1, size(is)
            i = is(k)
            z = B(k, l)

            found = .false.

            do m = A%g%ptr(j), A%g%ptr(j + 1) - 1
                if (A%g%node(m) == i) then
                    A%val(m) = z
                    found = .true.
                endif
            enddo

            if (.not. found) then
                call set_matrix_value_with_reallocation(A%g, A%val, j, i, z)
            endif
        enddo
   enddo

end subroutine csc_matrix_set_multiple_values



!--------------------------------------------------------------------------!
subroutine csc_matrix_add_multiple_values(A, is, js, B)                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csc_matrix), intent(inout) :: A
    integer, intent(in) :: is(:), js(:)
    real(dp), intent(in) :: B(:,:)
    ! local variables
    integer :: i, j, k, l, m
    real(dp) :: z
    logical :: found

    do l = 1, size(js)
        j = js(l)

        do k = 1, size(is)
            i = is(k)
            z = B(k, l)

            found = .false.

            do m = A%g%ptr(j), A%g%ptr(j + 1) - 1
                if (A%g%node(m) == i) then
                    A%val(m) = A%val(m) + z
                    found = .true.
                endif
            enddo

            if (.not. found) then
                call set_matrix_value_with_reallocation(A%g, A%val, j, i, z)
            endif
        enddo
   enddo

end subroutine csc_matrix_add_multiple_values



end module cs_matrices


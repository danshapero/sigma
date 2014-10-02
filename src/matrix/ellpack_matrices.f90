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
use graph_interfaces
use ellpack_graphs
use sparse_matrix_interfaces
use default_sparse_matrix_kernels,  only: get_degree_kernel, &
                        & get_degree_contiguous, get_degree_discontiguous

implicit none




!--------------------------------------------------------------------------!
type, extends(sparse_matrix_interface) :: ellpack_matrix                   !
!--------------------------------------------------------------------------!
    integer :: nnz
    class(ellpack_graph), pointer :: g => null()
    real(dp), allocatable :: val(:,:)
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: copy_graph => ellpack_matrix_copy_graph
    procedure :: set_graph => ellpack_matrix_set_graph


    !-----------
    ! Accessors
    !-----------
    procedure :: get_nnz => ellpack_matrix_get_nnz
    procedure :: get_value => ellpack_matrix_get_value
    procedure :: get_row_degree    => ellpack_matrix_get_row_degree
    procedure :: get_column_degree => ellpack_matrix_get_column_degree
    procedure :: get_row    => ellpack_matrix_get_row
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
    procedure :: set_dense_submatrix => ellpack_matrix_set_dense_submatrix
    procedure :: add_dense_submatrix => ellpack_matrix_add_dense_submatrix
    procedure :: zero => ellpack_matrix_zero
    procedure :: left_permute  => ellpack_matrix_left_permute
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


    !-------------------------
    ! Testing, debugging, I/O
    !-------------------------
    procedure :: to_dense_matrix => ellpack_matrix_to_dense_matrix


    !--------------------
    ! Auxiliary routines
    !--------------------
    procedure, private :: set_unallocated_matrix_value

end type ellpack_matrix




contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_matrix_copy_graph(A, g)                                 !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    class(graph_interface), intent(in) :: g

    if (A%nrow /= g%n .or. A%ncol /= g%m) then
        print *, 'Attemped to set ellpack matrix connectivity structure to'
        print *, 'graph with inconsistent dimensions.'
        print *, 'Dimensions of matrix:', A%nrow, A%ncol
        print *, 'Dimensions of graph: ', g%n, g%m
        call exit(1)
    endif

    if (.not. associated(A%g)) allocate(A%g)
    call A%g%copy(g)

    A%nnz = g%ne
    allocate(A%val(A%g%max_d, A%g%n))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine ellpack_matrix_copy_graph



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_set_graph(A, g)                                  !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    class(graph_interface), target, intent(in) :: g

    if (A%nrow /= g%n .or. A%ncol /= g%m) then
        print *, 'Attempted to set ellpack matrix connectivity structure to'
        print *, 'graph with inconsistent dimensions.'
        print *, 'Dimensions of matrix:', A%nrow, A%ncol
        print *, 'Dimensions of graph: ', g%n, g%m
        call exit(1)
    endif

    select type(g)
        class is(ellpack_graph)
            A%g => g
            call A%g%add_reference()
        class default
            print *, 'Attempted to set ellpack graph connectivity structure'
            print *, 'to point to a graph which is not an ellpack graph.'
            call exit(1)
    end select

    A%nnz = g%ne
    allocate(A%val(A%g%max_d, A%g%n))
    A%val = 0.0_dp

    A%graph_set = .true.

end subroutine ellpack_matrix_set_graph




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_matrix_get_nnz(A) result(nnz)                             !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    integer :: nnz

    nnz = A%nnz

end function ellpack_matrix_get_nnz



!--------------------------------------------------------------------------!
function ellpack_matrix_get_value(A, i, j) result(z)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in) :: i, j
    real(dp) :: z
    ! local variables
    integer :: k, d

    ! Set the return value to 0
    z = 0.0_dp

    d = A%g%degrees(i)
    do k = 1, d
        if (A%g%node(k, i) == j) z = A%val(k, i)
    enddo

end function ellpack_matrix_get_value



!--------------------------------------------------------------------------!
function ellpack_matrix_get_row_degree(A, k) result(d)                     !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in) :: k
    integer :: d

    d = A%g%degrees(k)

end function ellpack_matrix_get_row_degree



!--------------------------------------------------------------------------!
function ellpack_matrix_get_column_degree(A, k) result(d)                  !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in) :: k
    integer :: d

    d = get_degree_discontiguous(A%g, k)

end function ellpack_matrix_get_column_degree



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_get_row(A, nodes, slice, k)                      !
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

end subroutine ellpack_matrix_get_row



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_get_column(A, nodes, slice, k)                   !
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

end subroutine ellpack_matrix_get_column




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
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned

    call A%g%get_edges(edges, cursor, num_edges, num_returned)

end subroutine ellpack_matrix_get_edges



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_get_entries(A, edges, entries, cursor, &         !
                                                & num_edges, num_returned) !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    real(dp), intent(out) :: entries(num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned
    ! local variables
    integer :: i, i1, i2, num_added, num_from_this_row, idx

    associate(g => A%g)

    !--------
    ! Get a batch of values from `A` and store them in the array `entries`

    entries = 0.0_dp

    num_returned = min(num_edges, cursor%last - cursor%current)

    idx = cursor%idx
    i1 = cursor%current / g%max_d + 1
    i2 = (cursor%current + num_returned - 1) / g%max_d + 1

    num_added = 0

    do i = i1, i2
        num_from_this_row = min( g%max_d - idx, num_returned - num_added)

        entries(num_added + 1 : num_added + num_from_this_row) &
                        & = A%val( idx + 1 : idx + num_from_this_row, i)
        num_added = num_added + num_from_this_row

        idx = mod(idx + num_from_this_row, g%max_d)
    enddo


    !--------
    ! Get a batch of edges from the connectivity graph of `A`

    call g%get_edges(edges, cursor, num_edges, num_returned)

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
    integer :: k, d
    logical :: found

    found = .false.

    d = A%g%degrees(i)
    do k = 1, d
        if (A%g%node(k, i) == j) then
            A%val(k, i) = z
            found = .true.
        endif
    enddo

    if (.not. found) call A%set_unallocated_matrix_value(i, j, z)

end subroutine ellpack_matrix_set_value



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_add_value(A, i, j, z)                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: k, d
    logical :: found

    found = .false.

    d = A%g%degrees(i)
    do k = 1, d
        if (A%g%node(k, i) == j) then
            A%val(k, i) = A%val(k, i) + z
            found = .true.
        endif
    enddo

    if (.not. found) call A%set_unallocated_matrix_value(i, j, z)

end subroutine ellpack_matrix_add_value



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_set_dense_submatrix(A, is, js, B)                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: is(:), js(:)
    real(dp), intent(in) :: B(:,:)
    ! local variables
    integer :: i, j, k, l, m, d
    real(dp) :: z
    logical :: found

    do k = 1, size(is)
        i = is(k)
        d = A%g%degrees(i)

        do l = 1, size(js)
            j = js(l)
            z = B(k, l)

            found = .false.

            do m = 1, d
                if (A%g%node(m, i) == j) then
                    A%val(m, i) = z
                    found = .true.
                endif
            enddo

            if (.not. found) then
                call A%set_unallocated_matrix_value(i, j, z)
                d = A%g%degrees(i)
            endif
        enddo
    enddo

end subroutine ellpack_matrix_set_dense_submatrix



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_add_dense_submatrix(A, is, js, B)                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: is(:), js(:)
    real(dp), intent(in) :: B(:,:)
    ! local variables
    integer :: i, j, k, l, m, d
    real(dp) :: z
    logical :: found

    do k = 1, size(is)
        i = is(k)
        d = A%g%degrees(i)

        do l = 1, size(js)
            j = js(l)
            z = B(k, l)

            found = .false.

            do m = 1, d
                if (A%g%node(m, i) == j) then
                    A%val(m, i) = A%val(m, i) + z
                    found = .true.
                endif
            enddo

            if (.not. found) then
                call A%set_unallocated_matrix_value(i, j, z)
                d = A%g%degrees(i)
            endif
        enddo
    enddo

end subroutine ellpack_matrix_add_dense_submatrix



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_zero(A)                                          !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A

    A%val = 0.0_dp

end subroutine ellpack_matrix_zero



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_left_permute(A, p)                               !
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

end subroutine ellpack_matrix_left_permute



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_right_permute(A, p)                              !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine ellpack_matrix_right_permute




!==========================================================================!
!==== Matrix-vector multiplication                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_matvec_add(A, x, y)                                     !
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

end subroutine ellpack_matvec_add



!--------------------------------------------------------------------------!
subroutine ellpack_matvec_t_add(A, x, y)                                   !
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

end subroutine ellpack_matvec_t_add




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

end subroutine ellpack_matrix_destroy




!==========================================================================!
!==== Testing, debugging, I/O                                          ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_matrix_to_dense_matrix(A, B, trans)                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    real(dp), intent(out) :: B(:,:)
    logical, intent(in), optional :: trans
    ! local variables
    integer :: i, j, k, d, ind(2), ord(2)

    ! Set the dense matrix to 0
    B = 0.0_dp

    ! If we're actually making the transpose of A, set the variable `ord`
    ! so that we reverse the orientation of all edges
    ord = [1, 2]
    if (present(trans)) then
        if (trans) ord = [2, 1]
    endif

    do i = 1, A%nrow
        d = A%g%degrees(i)

        do k = 1, d
            j = A%g%node(k, i)

            ind = [i, j]
            ind = ind(ord)

            B(ind(1), ind(2)) = A%val(k, i)
        enddo
    enddo

end subroutine ellpack_matrix_to_dense_matrix




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

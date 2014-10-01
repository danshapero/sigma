!==========================================================================!
!==========================================================================!
module sparse_matrix_composites                                            !
!==========================================================================!
!==========================================================================!
!====     This module contains the sparse_matrix data type, which      ====!
!==== allows flexible wrapping and composition of the basic sparse     ====!
!==== matrix types.                                                    ====!
!====     The most basic usage of the sparse_matrix type is as a thin  ====!
!==== wrapper around the basic types. Rather than require the user to  ====!
!==== maintain a polymorphic                                           ====!
!====     class(sparse_matrix_interface), pointer :: A                 ====!
!==== object, she can instead interact with the non-polymorphic        ====!
!====     type(sparse_matrix) :: A                                     ====!
!==== object and request the desired storage format through methods    ====!
!==== rather than through typed allocation.                            ====!
!====     A sparse_matrix object can also be used as a composite of    ====!
!==== several sub-matrices, each of which could itself be a composite  ====!
!==== matrix. This functionality enables block decomposition of        ====!
!==== matrices, especially useful in multi-physics applications or     ====!
!==== higher-order finite element discretizations. It is also the      ====!
!==== primary means by which parallelism is enabled in SiGMA.          ====!
!==========================================================================!
!==========================================================================!

use graph_interfaces
use sparse_matrix_interfaces
use sparse_matrix_factory



!--------------------------------------------------------------------------!
type :: sparse_mat_pointer                                                 !
!--------------------------------------------------------------------------!
    class(sparse_matrix_interface), pointer :: mat => null()
end type sparse_mat_pointer



!--------------------------------------------------------------------------!
type, extends(sparse_matrix_interface) :: sparse_matrix                    !
!--------------------------------------------------------------------------!
    integer :: num_row_mats = 0, num_col_mats = 0
    integer, allocatable :: row_ptr(:), col_ptr(:)
    type(sparse_mat_pointer), allocatable, private :: sub_mats(:, :)

    logical, private :: num_blocks_set  = .false.
    logical, private :: block_sizes_set = .false.
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: set_dimensions => composite_mat_set_dimensions
    ! Set the global dimension of a composite matrix

    procedure :: set_num_blocks  => composite_mat_set_num_blocks
    ! Set the number of sub-matrices that the composite will be made of

    procedure :: set_block_sizes => composite_mat_set_block_sizes
    ! Set the size of each sub-matrix of the composite matrix

    procedure :: composite_mat_set_mat_type_leaf
    procedure :: composite_mat_set_mat_type_submat
    generic :: set_matrix_type => composite_mat_set_mat_type_leaf, &
                                & composite_mat_set_mat_type_submat
    ! Choose the sparse matrix format for each sub-matrix of the composite

    procedure :: copy_graph => composite_mat_copy_graph
    procedure :: set_graph  => composite_mat_set_graph
    ! Create the connectivity structure of a leaf matrix

    procedure :: copy_graph_submat => composite_mat_copy_graph_submat
    procedure :: set_graph_submat  => composite_mat_set_graph_submat
    ! Create the connectivity structure of a sub-matrix of the composite


    !-----------
    ! Accessors
    !-----------
    procedure :: get_nnz           => composite_mat_get_nnz
    procedure :: get_value         => composite_mat_get_value
    procedure :: get_submat_value  => composite_mat_get_submat_value
    procedure :: get_row_degree    => composite_mat_get_row_degree
    procedure :: get_column_degree => composite_mat_get_column_degree
    procedure :: get_row           => composite_mat_get_row
    procedure :: get_column        => composite_mat_get_column

    procedure :: is_leaf           => composite_mat_is_leaf
    ! Returns true if the sparse_matrix is a leaf, i.e. it wraps a single
    ! sparse matrix format and not a composite of several matrices

    procedure :: get_submatrix     => composite_mat_get_submatrix
    ! Return one of the composite's sub-matrices, wrapped in a sparse_matrix

    generic :: get => get_submat_value


    !-----------------------
    ! Edge, value iterators
    !-----------------------
    procedure :: make_cursor => composite_mat_make_cursor
    procedure :: get_edges   => composite_mat_get_edges
    procedure :: get_entries => composite_mat_get_entries


    !----------
    ! Mutators
    !----------
    procedure :: set_value           => composite_mat_set_value
    procedure :: add_value           => composite_mat_add_value
    ! procedure :: set_dense_submatrix => composite_mat_set_dense_submatrix
    ! procedure :: add_dense_submatrix => composite_mat_add_dense_submatrix
    procedure :: set_submat_value    => composite_mat_set_submat_value
    procedure :: add_submat_value    => composite_mat_add_submat_value
    procedure :: zero                => composite_mat_zero
    procedure :: left_permute        => composite_mat_left_permute
    procedure :: right_permute       => composite_mat_right_permute

    procedure :: set_submatrix       => composite_mat_set_submatrix

    ! Additional generics for mutators
    generic :: set => set_submat_value
    generic :: add => add_submat_value


    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure :: matvec_add   => composite_matvec_add
    procedure :: matvec_t_add => composite_matvec_t_add


    !-------------
    ! Destructors
    !-------------
    procedure :: destroy => composite_mat_destroy


    !----------------------
    ! Auxiliary procedures
    !----------------------
    procedure :: get_owning_row_matrix => composite_mat_get_owning_rowmat
    procedure :: get_owning_column_matrix => composite_mat_get_owning_colmat

end type sparse_matrix



private
public :: sparse_matrix



contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine composite_mat_set_dimensions(A, nrow, ncol)                     !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol

    A%nrow = nrow
    A%ncol = ncol

    A%dimensions_set = .true.

    call A%add_reference()

    if (A%is_leaf()) then
        call A%set_block_sizes([nrow], [ncol])
        call A%sub_mats(1, 1)%mat%set_dimensions(nrow, ncol)
    endif

end subroutine composite_mat_set_dimensions



!--------------------------------------------------------------------------!
subroutine composite_mat_set_num_blocks(A, num_row_mats, num_col_mats)     !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: num_row_mats, num_col_mats

    A%num_row_mats = num_row_mats
    A%num_col_mats = num_col_mats

    allocate( A%row_ptr(num_row_mats + 1) )
    allocate( A%col_ptr(num_col_mats + 1) )

    A%row_ptr = 0
    A%col_ptr = 0

    allocate(A%sub_mats(num_row_mats, num_col_mats))

    A%num_blocks_set = .true.

end subroutine composite_mat_set_num_blocks



!--------------------------------------------------------------------------!
subroutine composite_mat_set_block_sizes(A, rows, cols)                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:), cols(:)
    ! local variables
    integer :: it, jt

    if (.not. A%num_blocks_set) then
        call A%set_num_blocks(size(rows), size(cols))
    else
        if (size(rows) /= A%num_row_mats .or. &
            & size(cols) /= A%num_col_mats) then
            print *, "Attempted to set block sizes of a sparse matrix, "
            print *, "but arrays provided are of inconsistent dimensions"
            print *, "with already-set number of sub-matrices."
            print *, "Row, col array dimensions:", size(rows), size(cols)
            print *, "Number of row/col submatrices:", A%num_row_mats, &
                                                        & A%num_col_mats
            call exit(1)
        endif
    endif

    A%row_ptr(1) = 1
    A%col_ptr(1) = 1

    do it = 1, A%num_row_mats
        A%row_ptr(it + 1) = A%row_ptr(it) + rows(it)
    enddo

    do jt = 1, A%num_col_mats
        A%col_ptr(jt + 1) = A%col_ptr(jt) + cols(jt)
    enddo

    A%block_sizes_set = .true.

end subroutine composite_mat_set_block_sizes



!--------------------------------------------------------------------------!
subroutine composite_mat_set_mat_type_leaf(A, frmt)                        !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    character(len=*), intent(in) :: frmt

    if (.not. A%num_blocks_set) call A%set_num_blocks(1, 1)

    call choose_matrix_type(A%sub_mats(1, 1)%mat, frmt)

    if (A%dimensions_set) then
        call A%set_block_sizes([A%nrow], [A%ncol])
        call A%sub_mats(1, 1)%mat%set_dimensions(A%nrow, A%ncol)
    endif

end subroutine composite_mat_set_mat_type_leaf



!--------------------------------------------------------------------------!
subroutine composite_mat_set_mat_type_submat(A, it, jt, frmt)              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: it, jt
    character(len=*), intent(in) :: frmt

    if (.not. A%num_blocks_set) then
        print *, "Attempted to set the storage format of a sub-matrix"
        print *, "of a composite matrix which has not yet had the number"
        print *, "of sub-matrices set."
        call exit(1)
    endif

    call choose_matrix_type(A%sub_mats(it, jt)%mat, frmt)

    if (A%block_sizes_set) then
        call A%sub_mats(it, jt)%mat%set_dimensions( &
                & A%row_ptr(it + 1) - A%row_ptr(it), &
                & A%col_ptr(jt + 1) - A%col_ptr(jt))
    endif

end subroutine composite_mat_set_mat_type_submat



!--------------------------------------------------------------------------!
subroutine composite_mat_copy_graph(A, g)                                  !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    class(graph_interface), intent(in) :: g

    ! The matrix dimensions must have been set before giving it a graph
    ! structure
    if (.not. A%dimensions_set) then
        print *, "Attempted to copy a graph's structure to a composite "
        print *, "matrix which has not had its dimensions set."
        call exit(1)
    endif

    ! Ensure that we're not treating a composite matrix as if it were a leaf
    if (A%block_sizes_set .and. .not. A%is_leaf()) then
        print *, "Attempted to copy a graph's structure to a composite "
        print *, "matrix which is not a leaf."
        call exit(1)
    endif

    ! If the block structure has not already been set, then set it to
    ! reflect the fact that this matrix is a leaf
    if (.not. A%block_sizes_set) then
        call A%set_block_sizes( [A%nrow], [A%ncol] )
    endif

    ! Copy the structure of the graph `g` to the object that this matrix
    ! is wrapping
    call A%sub_mats(1, 1)%mat%copy_graph(g)

end subroutine composite_mat_copy_graph



!--------------------------------------------------------------------------!
subroutine composite_mat_set_graph(A, g)                                   !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    class(graph_interface), target, intent(in) :: g

    if (.not. A%dimensions_set) then
        print *, "Attempted to set the graph structure of a composite "
        print *, "matrix which has not had its dimensions set."
        call exit(1)
    endif

    if (A%block_sizes_set .and. .not. A%is_leaf()) then
        print *, "Attempted to set the graph structure of a composite "
        print *, "matrix which is not a leaf."
        call exit(1)
    endif

    if (.not. A%block_sizes_set) then
        call A%set_block_sizes( [A%nrow], [A%ncol] )
    endif

    call A%sub_mats(1, 1)%mat%set_graph(g)

end subroutine composite_mat_set_graph



!--------------------------------------------------------------------------!
subroutine composite_mat_copy_graph_submat(A, it, jt, g)                   !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: it, jt
    class(graph_interface), intent(in) :: g

    ! Check that the sub-matrix indices are in bounds
    if (it < 1 .or. it > A%num_row_mats &
        & .or. jt < 1 .or. jt > A%num_col_mats) then
        print *, "Sub-matrix indices", it, jt
        print *, "are out of bounds", A%num_row_mats, A%num_col_mats
        call exit(1)
    endif

    ! Copy the structure of `g` to the desired sub-matrix
    call A%sub_mats(it, jt)%mat%copy_graph(g)

end subroutine composite_mat_copy_graph_submat



!--------------------------------------------------------------------------!
subroutine composite_mat_set_graph_submat(A, it, jt, g)                    !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: it, jt
    class(graph_interface), intent(in) :: g

    if (it < 1 .or. it > A%num_row_mats &
        & .or. jt < 1 .or. jt > A%num_col_mats) then
        print *, "Sub-matrix indices", it, jt
        print *, "are out of bounds", A%num_row_mats, A%num_col_mats
        call exit(1)
    endif

    call A%sub_mats(it, jt)%mat%set_graph(g)

end subroutine composite_mat_set_graph_submat




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function composite_mat_get_nnz(A) result(nnz)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer :: nnz
    ! local variables
    integer :: it, jt

    nnz = 0

    do jt = 1, A%num_col_mats
        do it = 1, A%num_row_mats
            nnz = nnz + A%sub_mats(it, jt)%mat%get_nnz()
        enddo
    enddo

end function composite_mat_get_nnz


!--------------------------------------------------------------------------!
function composite_mat_get_value(A, i, j) result(z)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: i, j
    real(dp) :: z
    ! local variables
    integer :: ir, jr, it, jt
    class(sparse_matrix_interface), pointer :: C

    it = A%get_owning_row_matrix(i)
    jt = A%get_owning_column_matrix(j)

    C => A%sub_mats(it, jt)%mat

    ir = i - A%row_ptr(it) + 1
    jr = j - A%col_ptr(jt) + 1

    z = C%get_value(ir, jr)

end function composite_mat_get_value



!--------------------------------------------------------------------------!
function composite_mat_get_submat_value(A, it, jt, i, j) result(z)         !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: it, jt, i, j
    real(dp) :: z

    z = A%sub_mats(it, jt)%mat%get_value(i, j)

end function composite_mat_get_submat_value



!--------------------------------------------------------------------------!
function composite_mat_get_row_degree(A, k) result(d)                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: k
    integer :: d
    ! local variables
    integer :: ir, it, jt
    class(sparse_matrix_interface), pointer :: C

    it = A%get_owning_row_matrix(k)
    ir = k - A%row_ptr(it) + 1

    d = 0
    do jt = 1, A%num_col_mats
        C => A%sub_mats(it, jt)%mat
        d = d + C%get_row_degree(ir)
    enddo

end function composite_mat_get_row_degree



!--------------------------------------------------------------------------!
function composite_mat_get_column_degree(A, k) result(d)                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: k
    integer :: d
    ! local variables
    integer :: it, jt, jr
    class(sparse_matrix_interface), pointer :: C

    jt = A%get_owning_column_matrix(k)
    jr = k - A%col_ptr(jt) + 1

    d = 0
    do it = 1, A%num_row_mats
        C => A%sub_mats(it, jt)%mat
        d = d + C%get_column_degree(jr)
    enddo

end function composite_mat_get_column_degree



!--------------------------------------------------------------------------!
subroutine composite_mat_get_row(A, nodes, slice, k)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k
    ! local variables
    integer :: ir, it, jt, l, d
    class(sparse_matrix_interface), pointer :: C

    it = A%get_owning_row_matrix(k)
    ir = k - A%row_ptr(it) + 1

    l = 0
    d = 0

    do jt = 1, A%num_col_mats
        C => A%sub_mats(it, jt)%mat

        d = C%get_row_degree(ir)
        call C%get_row(nodes(l + 1 : l + d), slice(l + 1 : l + d), ir)

        nodes(l + 1 : l + d) = nodes(l + 1 : l + d) + A%col_ptr(jt) - 1

        l = l + d
    enddo

end subroutine composite_mat_get_row



!--------------------------------------------------------------------------!
subroutine composite_mat_get_column(A, nodes, slice, k)                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k
    ! local variables
    integer :: it, jt, jr, l, d
    class(sparse_matrix_interface), pointer :: C

    jt = A%get_owning_column_matrix(k)
    jr = k - A%col_ptr(jt) + 1

    l = 0
    d = 0

    do it = 1, A%num_row_mats
        C => A%sub_mats(it, jt)%mat

        d = C%get_column_degree(jr)
        call C%get_column(nodes(l+1 : l+d), slice(l+1 : l+d), jr)

        nodes(l+1 : l+d) = nodes(l+1 : l+d) + A%row_ptr(it) - 1

        l = l + d
    enddo

end subroutine composite_mat_get_column



!--------------------------------------------------------------------------!
function composite_mat_is_leaf(A) result(leaf)                             !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    logical :: leaf

    leaf = A%num_row_mats == 1 .and. A%num_col_mats == 1

end function composite_mat_is_leaf



!--------------------------------------------------------------------------!
function composite_mat_get_submatrix(A, it, jt) result(C)                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: it, jt
    type(sparse_matrix), pointer :: C
    ! local variables
    class(sparse_matrix_interface), pointer :: Aij

    Aij => A%sub_mats(it, jt)%mat

    ! Check whether the desired sub-matrix is a leaf or a composite matrix.
    select type(Aij)
        ! If it's a composite matrix, we can return a pointer to it
        type is(sparse_matrix)
            C => Aij

        ! Otherwise, make C a leaf matrix referring to the desired sub-
        ! matrix and return a pointer to C.
        class default
            allocate(C)
            call C%set_dimensions(A%row_ptr(it+1) - A%row_ptr(it), &
                                & A%col_ptr(jt+1) - A%col_ptr(jt))

            call C%set_block_sizes([C%nrow], [C%ncol])

            C%sub_mats(1, 1)%mat => Aij

    end select

    call Aij%add_reference()

end function composite_mat_get_submatrix




!==========================================================================!
!==== Edge, value iterators                                            ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function composite_mat_make_cursor(A) result(cursor)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    type(graph_edge_cursor) :: cursor
    ! local variables
    integer :: it, jt

    cursor%first = 1
    cursor%last = A%get_nnz()
    cursor%current = 0
    cursor%edge = [1, 1]

    allocate(cursor%sub_cursor(A%num_row_mats, A%num_col_mats))

    do jt = 1, A%num_col_mats
        do it = 1, A%num_row_mats
            cursor%sub_cursor(it, jt) = A%sub_mats(it, jt)%mat%make_cursor()
        enddo
    enddo

end function composite_mat_make_cursor



!--------------------------------------------------------------------------!
subroutine composite_mat_get_edges(A, edges, cursor, &                     !
                                                & num_edges, num_returned) !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned
    ! local variables
    integer :: it, jt
    class(sparse_matrix_interface), pointer :: C

    it = cursor%edge(1)
    jt = cursor%edge(2)

    C => A%sub_mats(it, jt)%mat

    call C%get_edges(edges, cursor%sub_cursor(it, jt), &
                                                  & num_edges, num_returned)

    cursor%current = cursor%current + num_returned

    ! Some logic to increment it/jt if we've reached the end of that matrix
    if (cursor%sub_cursor(it, jt)%current == &
                                    & cursor%sub_cursor(it, jt)%last) then
        it = mod(it, A%num_row_mats) + 1
        if (it == A%num_row_mats) jt = jt + 1
    endif

    cursor%edge(1) = it
    cursor%edge(2) = jt

    if (cursor%current == cursor%last) deallocate(cursor%sub_cursor)

end subroutine composite_mat_get_edges



!--------------------------------------------------------------------------!
subroutine composite_mat_get_entries(A, edges, entries, cursor, &          !
                                                & num_edges, num_returned) !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    real(dp), intent(out) :: entries(num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned
    ! local variables
    integer :: it, jt
    class(sparse_matrix_interface), pointer :: C

    it = cursor%edge(1)
    jt = cursor%edge(2)

    C => A%sub_mats(it, jt)%mat

    call C%get_entries(edges, entries, cursor%sub_cursor(it, jt), &
                                                & num_edges, num_returned)

    cursor%current = cursor%current + num_returned

    if (cursor%sub_cursor(it, jt)%current == &
                                    & cursor%sub_cursor(it, jt)%last) then
        it = mod(it, A%num_row_mats) + 1
        if (it == A%num_row_mats) jt = jt + 1
    endif

    cursor%edge(1) = it
    cursor%edge(2) = jt

end subroutine composite_mat_get_entries




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine composite_mat_set_value(A, i, j, z)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: ir, jr, it, jt
    class(sparse_matrix_interface), pointer :: C

    it = A%get_owning_row_matrix(i)
    jt = A%get_owning_column_matrix(j)

    ir = i - A%row_ptr(it) + 1
    jr = j - A%col_ptr(jt) + 1

    C => A%sub_mats(it, jt)%mat

    call C%set_value(ir, jr, z)

end subroutine composite_mat_set_value



!--------------------------------------------------------------------------!
subroutine composite_mat_add_value(A, i, j, z)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    integer :: ir, jr, it, jt
    class(sparse_matrix_interface), pointer :: C

    it = A%get_owning_row_matrix(i)
    jt = A%get_owning_column_matrix(j)

    ir = i - A%row_ptr(it) + 1
    jr = j - A%col_ptr(jt) + 1

    C => A%sub_mats(it, jt)%mat

    call C%add_value(ir, jr, z)

end subroutine composite_mat_add_value



!--------------------------------------------------------------------------!
subroutine composite_mat_set_submat_value(A, it, jt, i, j, z)              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: it, jt, i, j
    real(dp), intent(in) :: z

    call A%sub_mats(it, jt)%mat%set_value(i, j, z)

end subroutine composite_mat_set_submat_value



!--------------------------------------------------------------------------!
subroutine composite_mat_add_submat_value(A, it, jt, i, j, z)              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: it, jt, i, j
    real(dp), intent(in) :: z

    call A%sub_mats(it, jt)%mat%add_value(i, j, z)

end subroutine composite_mat_add_submat_value



!--------------------------------------------------------------------------!
subroutine composite_mat_zero(A)                                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    ! local variables
    integer :: it, jt
    class(sparse_matrix_interface), pointer :: C

    do jt = 1, A%num_col_mats
        do it = 1, A%num_row_mats
            C => A%sub_mats(it, jt)%mat
            call C%zero()
        enddo
    enddo

end subroutine composite_mat_zero



!--------------------------------------------------------------------------!
subroutine composite_mat_left_permute(A, p)                                !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

	if (A%is_leaf()) then
		call A%sub_mats(1, 1)%mat%left_permute(p)
	else
		print *, 'Cannot permute a composite matrix.'
		call exit(1)
	endif

end subroutine composite_mat_left_permute



!--------------------------------------------------------------------------!
subroutine composite_mat_right_permute(A, p)                               !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

	if (A%is_leaf()) then
		call A%sub_mats(1, 1)%mat%right_permute(p)
	else
		print *, 'Cannot permute a composite matrix.'
		call exit(1)
	endif

end subroutine composite_mat_right_permute



!--------------------------------------------------------------------------!
subroutine composite_mat_set_submatrix(A, it, jt, B)                       !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: it, jt
    class(sparse_matrix_interface), target :: B

    select type(B)
        class is(sparse_matrix)
            if (B%is_leaf()) then
                A%sub_mats(it, jt)%mat => B%sub_mats(1, 1)%mat
            endif
        class default
            A%sub_mats(it, jt)%mat => B
    end select

    call A%sub_mats(it, jt)%mat%add_reference()

end subroutine composite_mat_set_submatrix




!==========================================================================!
!==== Matrix-vector multiplication                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine composite_matvec_add(A, x, y)                                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i1, i2, j1, j2, it, jt
    class(sparse_matrix_interface), pointer :: C

    do it = 1, A%num_row_mats     ! This loop can be parallelized
        i1 = A%row_ptr(it)
        i2 = A%row_ptr(it + 1) - 1

        do jt = 1, A%num_col_mats
            j1 = A%col_ptr(jt)
            j2 = A%col_ptr(jt + 1) - 1

            C => A%sub_mats(it, jt)%mat

            call C%matvec_add( x(j1 : j2), y(i1 : i2) )
        enddo
    enddo

end subroutine composite_matvec_add



!--------------------------------------------------------------------------!
subroutine composite_matvec_t_add(A, x, y)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i1, i2, j1, j2, it, jt
    class(sparse_matrix_interface), pointer :: C

    do jt = 1, A%num_col_mats     ! This loop can be parallelized
        j1 = A%col_ptr(jt)
        j2 = A%col_ptr(jt + 1) - 1

        do it = 1, A%num_row_mats
            i1 = A%row_ptr(it)
            i2 = A%row_ptr(it + 1) - 1

            C => A%sub_mats(it, jt)%mat

            call C%matvec_t_add( x(i1 : i2), y(j1 : j2) )
        enddo
    enddo

end subroutine composite_matvec_t_add




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine composite_mat_destroy(A)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    ! local variables
    integer :: it, jt
    class(sparse_matrix_interface), pointer :: C

    do jt = 1, A%num_col_mats
        do it = 1, A%num_row_mats
            C => A%sub_mats(it, jt)%mat

            if (associated(C)) then
                call C%remove_reference()

                if (C%reference_count <= 0) then
                    call C%destroy()
                    deallocate(C)
                endif
            endif

            nullify(A%sub_mats(it, jt)%mat)
        enddo
    enddo

    deallocate(A%row_ptr)
    deallocate(A%col_ptr)
    deallocate(A%sub_mats)
    A%num_row_mats = 0
    A%num_col_mats = 0

    A%reference_count = 0

    A%nrow = 0
    A%ncol = 0

    A%dimensions_set = .false.
    A%graph_set = .false.
    A%num_blocks_set = .false.
    A%block_sizes_set = .false.

end subroutine composite_mat_destroy



!==========================================================================!
!==== Auxiliary procedures                                             ====!
!==========================================================================!

!--------------------------------------------------------------------------!
elemental function composite_mat_get_owning_rowmat(A, i) result(it)        !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: i
    integer :: it

    do it = 1, A%num_row_mats
        if (A%row_ptr(it) <= i .and. A%row_ptr(it + 1) > i) exit
    enddo

end function composite_mat_get_owning_rowmat



!--------------------------------------------------------------------------!
elemental function composite_mat_get_owning_colmat(A, j) result(jt)        !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: j
    integer :: jt

    do jt = 1, A%num_col_mats
        if (A%col_ptr(jt) <= j .and. A%col_ptr(jt + 1) > j) exit
    enddo

end function composite_mat_get_owning_colmat




end module sparse_matrix_composites


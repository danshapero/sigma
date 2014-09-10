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
    integer :: num_row_mats, num_col_mats
    integer, allocatable :: row_ptr(:), col_ptr(:)
    type(sparse_mat_pointer), allocatable, private :: sub_mats(:, :)
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: copy_graph_structure => composite_mat_copy_graph_structure
    procedure :: set_graph            => composite_mat_set_graph

    procedure :: copy_graph_submat => composite_mat_copy_graph_submat
    procedure :: set_graph_submat  => composite_mat_set_graph_submat
    procedure :: set_block_structure => composite_mat_set_block_structure


    !-----------
    ! Accessors
    !-----------
    procedure :: get_value         => composite_mat_get_value
    procedure :: get_row_degree    => composite_mat_get_row_degree
    procedure :: get_column_degree => composite_mat_get_column_degree
    procedure :: get_row           => composite_mat_get_row
    procedure :: get_column        => composite_mat_get_column

    procedure :: is_leaf           => composite_mat_is_leaf
    procedure :: get_submatrix     => composite_mat_get_submatrix


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
    procedure :: zero                => composite_mat_zero
    procedure :: left_permute        => composite_mat_left_permute
    procedure :: right_permute       => composite_mat_right_permute

    procedure :: set_submatrix       => composite_mat_set_submatrix


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
subroutine composite_mat_copy_graph_structure(A, g)                        !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    class(graph_interface), intent(in) :: g

    if (.not. A%is_leaf()) call exit(1)

    call A%sub_mats(1, 1)%mat%copy_graph_structure(g)

end subroutine composite_mat_copy_graph_structure



!--------------------------------------------------------------------------!
subroutine composite_mat_set_graph(A, g)                                   !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    class(graph_interface), target, intent(in) :: g

    if (.not. A%is_leaf()) call exit(1)

    call A%sub_mats(1, 1)%mat%set_graph(g)

end subroutine composite_mat_set_graph



!--------------------------------------------------------------------------!
subroutine composite_mat_copy_graph_submat(A, it, jt, g)                   !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: it, jt
    class(graph_interface), intent(in) :: g

    call A%sub_mats(it, jt)%mat%copy_graph_structure(g)

end subroutine composite_mat_copy_graph_submat



!--------------------------------------------------------------------------!
subroutine composite_mat_set_graph_submat(A, it, jt, g)                    !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: it, jt
    class(graph_interface), intent(in) :: g

    call A%sub_mats(it, jt)%mat%set_graph(g)

end subroutine composite_mat_set_graph_submat



!--------------------------------------------------------------------------!
subroutine composite_mat_set_block_structure(A, rows, cols)                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:), cols(:)
    ! local variables
    integer :: it, jt

    A%num_row_mats = size(rows)
    A%num_col_mats = size(cols)

    allocate( A%row_ptr(A%num_row_mats + 1) )
    allocate( A%col_ptr(A%num_col_mats + 1) )

    A%row_ptr(1) = 1
    A%col_ptr(1) = 1

    do it = 1, A%num_row_mats
        A%row_ptr(it + 1) = A%row_ptr(it) + rows(it)
    enddo

    do jt = 1, A%num_col_mats
        A%col_ptr(jt + 1) = A%col_ptr(jt) + cols(jt)
    enddo

    allocate(A%sub_mats(it, jt))

end subroutine composite_mat_set_block_structure




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

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

            C%num_row_mats = 1
            C%num_col_mats = 1

            allocate(C%row_ptr(2), C%col_ptr(2))

            C%row_ptr(1) = 1
            C%row_ptr(2) = C%nrow + 1

            C%col_ptr(1) = 1
            C%col_ptr(2) = C%ncol + 1

            allocate(C%sub_mats(1, 1))
            C%sub_mats(1, 1)%mat => Aij

            call Aij%add_reference()
    end select

end function composite_mat_get_submatrix




!==========================================================================!
!==== Edge, value iterators                                            ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function composite_mat_make_cursor(A) result(cursor)                       !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    type(graph_edge_cursor) :: cursor


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

    do it = 1, A%num_row_mats
        do jt = 1, A%num_col_mats
            C => A%sub_mats(it, jt)%mat

            call C%remove_reference()

            if (C%reference_count <= 0) then
                call C%destroy()
                deallocate(C)
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
    A%nnz = 0

    A%dimensions_set = .false.
    A%graph_set = .false.

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


!==========================================================================!
!==========================================================================!
module sparse_matrix_composites                                            !
!==========================================================================!
!==========================================================================!

use sparse_matrix_interfaces

implicit none



!--------------------------------------------------------------------------!
type :: sparse_mat_pointer                                                 !
!--------------------------------------------------------------------------!
    class(sparse_matrix_interface), pointer :: mat
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

    procedure :: set_block_structure => composite_mat_set_block_structure
    procedure :: set_storage_format =>  composite_mat_set_storage_format


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
    procedure :: set_dense_submatrix => composite_mat_set_dense_submatrix
    procedure :: add_dense_submatrix => composite_mat_add_dense_submatrix
    procedure :: zero                => composite_mat_zero
    procedure :: left_permute        => composite_mat_left_permute
    procedure :: right_permute       => composite_mat_right_permute

    procedure :: set_submatrix


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



contains



!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!



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

    it = A%get_owning_rowmat(i)
    jt = A%get_owning_colmat(j)

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
    ir = i - A%row_ptr(it) + 1

    d = 0
    do jt = 1, num_col_mats
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
    jr = j - A%col_ptr(jt) + 1

    d = 0
    do it = 1, num_row_mats
        C => A%sub_mats(it, jt)%mat
        d = d + C%get_column_degree(jr)
    enddo

end function composite_mat_get_column_degree




!==========================================================================!
!==== Edge, value iterators                                            ====!
!==========================================================================!




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine composite_mat_set_value(A, i, j, z)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer :: i, j
    real(dp) :: z
    ! local variables
    integer :: ir, jr, it, jt

end subroutine composite_mat_set_value



!--------------------------------------------------------------------------!
subroutine composite_mat_add_value(A, i, j, z)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer :: i, j
    real(dp) :: z
    ! local variables
    integer :: ir, jr, it, jt

end subroutine composite_mat_add_value




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

    do it = 1, num_row_mats     ! This loop can be parallelized
        i1 = row_ptr(it)
        i2 = row_ptr(it + 1) - 1

        do jt = 1, num_col_mats
            j1 = col_ptr(jt)
            j2 = col_ptr(jt + 1) - 1

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

    do jt = 1, num_col_mats     ! This loop can be parallelized
        j1 = col_ptr(jt)
        j2 = col_ptr(jt + 1) - 1

        do it = 1, num_row_mats
            i1 = row_ptr(it)
            i2 = row_ptr(it + 1) - 1

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

    do it = 1, num_row_mats
        do jt = 1, num_col_mats
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
elemental function sparse_mat_get_owning_rowmat(A, i) result(it)           !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: i
    integer, intent(out) :: it

    do it = 1, num_row_mats
        if (A%row_ptr(it) <= i .and. A%row_ptr(it + 1) > i) exit
    enddo

end function sparse_mat_get_owning_rowmat



!--------------------------------------------------------------------------!
elemental function sparse_mat_get_owning_colmat(A, j) result(jt)           !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: j
    integer, intent(out) :: jt

    do jt = 1, num_col_mats
        if (A%col_ptr(jt) <= j .and. A%col_ptr(jt + 1) > j) exit
    enddo

end function sparse_mat_get_owning_colmat




end module sparse_matrix_composites


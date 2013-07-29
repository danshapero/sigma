module sparse_matrices

use graphs

implicit none



!--------------------------------------------------------------------------!
type :: sparse_matrix                                                      !
!--------------------------------------------------------------------------!
    ! Data for sparse matrix storage
    type(graph), pointer :: g
    real(kind(1d0)), allocatable :: val(:)
    logical :: symmetric, positive
    ! Object-bound procedures for matrix operations
    procedure(get_value_ifc), pointer       :: get_value
    procedure(set_value_ifc), pointer       :: set_value
end type sparse_matrix



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    function get_value_ifc(A,i,j) result(val)
        import :: sparse_matrix
        class(sparse_matrix), intent(in) :: A
        integer, intent(in) :: i,j
        real(kind(1d0)) :: val
    end function get_value_ifc

    subroutine set_value_ifc(A,i,j,val)
        import :: sparse_matrix
        class(sparse_matrix), intent(inout) :: A
        integer, intent(in) :: i,j
        real(kind(1d0)), intent(in) :: val
    end subroutine set_value_ifc

end interface





contains



function get_sparse_matrix(storage,g) result(A)
    implicit none
    character(len=*), intent(in), optional :: storage
    type(graph), pointer, intent(in), optional :: g
    type(sparse_matrix) :: A

    if (present(g)) then
        A%g => g
    else
        allocate(A%g)
        g = get_graph(storage)
    endif

end function get_sparse_matrix



!--------------------------------------------------------------------------!
! Sparse matrix operations for the coordinate format                       !
!--------------------------------------------------------------------------!

function coo_get_value(A,i,j) result(val)
    implicit none
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)) :: val
    ! local variables
    integer :: n

    val = 0.d0

    do n=1,A%g%ne
        if ( A%g%ia(n)==i .and. A%g%ja(n)==j ) val = A%val(A%g%nrow+n)
    enddo

end function coo_get_value



subroutine coo_set_value(A,i,j,val)
    implicit none
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: val
    ! local variables
    integer :: n

    do n=1,A%g%ne
        if ( A%g%ia(n)==i .and. A%g%ja(n)==j ) A%val(A%g%nrow+n) = val
    enddo

end subroutine coo_set_value



end module sparse_matrices

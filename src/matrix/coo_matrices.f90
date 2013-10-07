module coo_matrices

use sparse_matrices
use coo_graphs

implicit none



!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: coo_matrix                                 !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    class(coo_graph), pointer :: g
contains
    procedure :: init => coo_matrix_init
    procedure :: neighbors => coo_matrix_neighbors
    procedure :: get_value => coo_get_value
    procedure :: set_value => coo_set_value
    procedure :: add_value => coo_add_value
    procedure :: sub_matrix_add => coo_sub_matrix_add
    procedure :: left_permute => coo_left_permute
    procedure :: right_permute => coo_right_permute
    procedure :: matvec => coo_matvec
    procedure :: matvec_t => coo_matvec_t
    procedure :: matvec_add => coo_matvec_add
    procedure :: matvec_t_add => coo_matvec_t_add
    procedure, private :: coo_set_value_not_preallocated
end type coo_matrix




contains





!--------------------------------------------------------------------------!
subroutine coo_matrix_init(A,nrow,ncol,orientation,g)                      !
!--------------------------------------------------------------------------!
    class(coo_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol
    character(len=3), intent(in) :: orientation
    class(graph), pointer, intent(in), optional :: g

    A%nrow = nrow
    A%ncol = ncol
    A%orientation = orientation

    if (present(g)) then
        select type(g)
            class is(coo_graph)
                A%g => g
            class default
                print *, 'Structure graph g of COO matrix A must be '
                print *, 'a COO graph. Exiting.'
                call exit(1)
        end select

        A%nrow = g%n
        A%ncol = g%n
    else
        allocate(coo_graph::A%g)
        call A%g%init(nrow,ncol)
    endif

    A%nnz = A%g%ne
    allocate(A%val(A%g%edges(1)%capacity))
    A%max_degree = A%g%max_degree

end subroutine coo_matrix_init



!--------------------------------------------------------------------------!
subroutine coo_matrix_neighbors(A,i,nbrs)                                  !
!--------------------------------------------------------------------------!
    class(coo_matrix), intent(in) :: A
    integer, intent(in)  :: i
    integer, intent(out) :: nbrs(:)

    nbrs = 0
    call A%g%neighbors(i,nbrs)

end subroutine coo_matrix_neighbors



!--------------------------------------------------------------------------!
function coo_get_value(A,i,j)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: coo_get_value
    ! local variables
    integer :: k

    coo_get_value = 0.0_dp
    k = A%g%find_edge(i,j)
    if (k/=-1) coo_get_value = A%val(k)

end function coo_get_value



!--------------------------------------------------------------------------!
subroutine coo_set_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%g%find_edge(i,j)
    if (k/=-1) then
        A%val(k) = val
    else
        call A%coo_set_value_not_preallocated(i,j,val)
    endif

end subroutine coo_set_value



!--------------------------------------------------------------------------!
subroutine coo_add_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%g%find_edge(i,j)
    if (k/=-1) then
        A%val(k) = A%val(k)+val
    else
        call A%coo_set_value_not_preallocated(i,j,val)
    endif

end subroutine coo_add_value



!--------------------------------------------------------------------------!
subroutine coo_sub_matrix_add(A,B)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_matrix), intent(inout) :: A
    class(sparse_matrix), intent(in) :: B
    ! local variables
    integer :: i,j,k

    do k=1,A%nnz
        i = A%g%edges(1)%get_entry(k)
        j = A%g%edges(2)%get_entry(k)
        A%val(k) = A%val(k)+B%get_value(i,j)
    enddo

end subroutine coo_sub_matrix_add



!--------------------------------------------------------------------------!
subroutine coo_left_permute(A,p)                                           !
!--------------------------------------------------------------------------!
    class(coo_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%left_permute(p)

end subroutine coo_left_permute



!--------------------------------------------------------------------------!
subroutine coo_right_permute(A,p)                                          !
!--------------------------------------------------------------------------!
    class(coo_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine coo_right_permute



!--------------------------------------------------------------------------!
subroutine coo_matvec_add(A,x,y)                                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i,j,k

    do k=1,A%nnz
        i = A%g%edges(1)%get_entry(k)
        j = A%g%edges(2)%get_entry(k)
        y(i) = y(i)+A%val(k)*x(j)
    enddo

end subroutine coo_matvec_add



!--------------------------------------------------------------------------!
subroutine coo_matvec_t_add(A,x,y)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i,j,k

    do k=1,A%nnz
        i = A%g%edges(2)%get_entry(k)
        j = A%g%edges(1)%get_entry(k)
        y(i) = y(i)+A%val(k)*x(j)
    enddo

end subroutine coo_matvec_t_add



!--------------------------------------------------------------------------!
subroutine coo_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    class(coo_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    y = 0.0_dp
    call A%matvec_add(x,y)

end subroutine coo_matvec



!--------------------------------------------------------------------------!
subroutine coo_matvec_t(A,x,y)                                             !
!--------------------------------------------------------------------------!
    class(coo_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    y = 0.0_dp
    call A%matvec_t_add(x,y)

end subroutine coo_matvec_t



!--------------------------------------------------------------------------!
subroutine coo_set_value_not_preallocated(A,i,j,val)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: capacity
    real(dp), allocatable :: val_temp(:)

    capacity = A%g%edges(1)%capacity

    ! If A is at capacity, we need to increase the storage size
    if (A%nnz==capacity) then
        ! Make a temporary array with twice the storage capacity as the
        ! array for the values of A
        allocate(val_temp(2*capacity))

        ! Copy the values of A over into the temporary array
        val_temp(1:A%nnz) = A%val(1:A%nnz)

        ! Move the allocation from the temporary array to the values of A
        call move_alloc(from=val_temp, to=A%val)
    endif

    ! Add the new edge to the graph g storing A's structure
    call A%g%add_edge(i,j)

    ! Add in the new matrix entry
    A%val(A%nnz+1) = val

    ! Increment the number of non-zero entries
    A%nnz = A%nnz+1

    ! Increment the maximum degree of the matrix, if need be
    A%max_degree = A%g%max_degree

end subroutine coo_set_value_not_preallocated



end module coo_matrices

module csr_matrices

use sparse_matrices
use cs_graphs

implicit none



!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: csr_matrix                                 !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    class(cs_graph), pointer :: g
contains
    procedure :: init => csr_init
    procedure :: assemble => csr_assemble
    procedure :: neighbors => csr_neighbors
    procedure :: get_value => csr_get_value
    procedure :: set_value => csr_set_value
    procedure :: add_value => csr_add_value
    procedure :: matvec => csr_matvec
    procedure, private :: csr_set_value_not_preallocated
end type csr_matrix




contains





!--------------------------------------------------------------------------!
subroutine csr_init(A,nrow,ncol)                                           !
!--------------------------------------------------------------------------!
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol

    A%nrow = nrow
    A%ncol = ncol
    A%max_degree = 0

end subroutine csr_init



!--------------------------------------------------------------------------!
subroutine csr_assemble(A,g)                                               !
!--------------------------------------------------------------------------!
    class(csr_matrix), intent(inout) :: A
    class(cs_graph), pointer, intent(in) :: g

    A%g => g

    A%nrow = g%n
    A%ncol = g%m
    A%nnz = g%ne
    A%max_degree = g%max_degree

    allocate(A%val(A%nnz))

end subroutine csr_assemble



!--------------------------------------------------------------------------!
subroutine csr_neighbors(A,i,nbrs)                                         !
!--------------------------------------------------------------------------!
    class(csr_matrix), intent(in) :: A
    integer, intent(in)  :: i
    integer, intent(out) :: nbrs(:)

    nbrs = 0
    call A%g%neighbors(i,nbrs)

end subroutine csr_neighbors



!--------------------------------------------------------------------------!
function csr_get_value(A,i,j)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: csr_get_value
    ! local variables
    integer :: k

    csr_get_value = 0_dp
    k = A%g%find_edge(i,j)
    if (k/=-1) csr_get_value = A%val(k)

end function csr_get_value



!--------------------------------------------------------------------------!
subroutine csr_set_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%g%find_edge(i,j)
    if (k/=-1) then
        A%val(k) = val
    else
        call A%csr_set_value_not_preallocated(i,j,val)
    endif

end subroutine csr_set_value



!--------------------------------------------------------------------------!
subroutine csr_add_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%g%find_edge(i,j)
    if (k/=-1) then
        A%val(k) = A%val(k)+val
    else
        call A%csr_set_value_not_preallocated(i,j,val)
    endif

end subroutine csr_add_value



!--------------------------------------------------------------------------!
subroutine csr_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z

    do i=1,A%nrow
        z = 0_dp
        do k=A%g%ia(i),A%g%ia(i+1)-1
            j = A%g%ja(k)
            z = z+A%val(k)*x(j)
        enddo
        y(i) = z
    enddo

end subroutine csr_matvec



!--------------------------------------------------------------------------!
subroutine csr_set_value_not_preallocated(A,i,j,val)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    real(dp) :: val_temp(A%nnz)
    integer :: k

    call A%g%add_edge(i,j)
    k = A%g%find_edge(i,j)
    val_temp = A%val
    deallocate(A%val)
    allocate(A%val(A%nnz+1))
    A%val(1:k-1) = val_temp(1:k-1)
    A%val(k) = val
    A%val(k+1:A%nnz+1) = val_temp(k:A%nnz)
    A%nnz = A%nnz+1
    A%max_degree = A%g%max_degree

end subroutine csr_set_value_not_preallocated




end module csr_matrices

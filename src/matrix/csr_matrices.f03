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
    procedure :: get_value => csr_get_value
    procedure :: set_value => csr_set_value
    procedure :: add_value => csr_add_value
    procedure :: matvec => csr_matvec
end type csr_matrix




contains





!--------------------------------------------------------------------------!
subroutine csr_init(A,nrow,ncol)                                           !
!--------------------------------------------------------------------------!
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol

    A%nrow = nrow
    A%ncol = ncol

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

    allocate(A%val(A%nnz))

end subroutine csr_assemble



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
    if (A%g%connected(i,j)) then
        k = A%g%find_edge(i,j)
        csr_get_value = A%val(k)
    endif

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

    if (.not.A%g%connected(i,j)) then
        call A%g%add_edge(i,j)
    endif

    k = A%g%find_edge(i,j)
    A%val(k) = val

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

    if (.not.A%g%connected(i,j)) then
        call A%g%add_edge(i,j)
    endif

    k = A%g%find_edge(i,j)
    A%val(k) = A%val(k)+val

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






end module csr_matrices

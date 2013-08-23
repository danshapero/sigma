module cs_matrices

use sparse_matrices
use cs_graphs

implicit none



!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: cs_matrix                                  !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    class(cs_graph), pointer :: g
    ! procedure pointers to implementations of matrix operations
    procedure(cs_find_entry_ifc), pointer :: find_entry
    procedure(cs_permute_ifc), pointer :: left_permute_impl, &
                                        & right_permute_impl
    procedure(cs_matvec_ifc), pointer :: matvec_impl, matvec_t_impl
contains
    ! front-ends to matrix operations
    procedure :: init => cs_matrix_init
    procedure :: assemble => cs_assemble
    procedure :: neighbors => cs_matrix_neighbors
    procedure :: get_value => cs_get_value
    procedure :: set_value => cs_set_value, add_value => cs_add_value
    procedure :: sub_matrix_add => cs_sub_matrix_add
    procedure :: left_permute => cs_left_permute, &
                & right_permute => cs_right_permute
    procedure :: matvec => cs_matvec, matvec_t => cs_matvec_t
    procedure, private :: cs_set_value_not_preallocated
end type cs_matrix




!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    function cs_find_entry_ifc(A,i,j)
        import :: cs_matrix
        class(cs_matrix), intent(in) :: A
        integer, intent(in) :: i,j
        integer :: cs_find_entry_ifc
    end function cs_find_entry_ifc

    subroutine cs_matvec_ifc(A,x,y)
        import :: cs_matrix, dp
        class(cs_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
    end subroutine

    subroutine cs_permute_ifc(A,p)
        import :: cs_matrix
        class(cs_matrix), intent(inout) :: A
        integer, intent(in) :: p(:)
    end subroutine cs_permute_ifc

end interface




contains



!==========================================================================!
!==========================================================================!
!== Front-ends to matrix operations                                      ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine cs_matrix_init(A,nrow,ncol)                                     !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol

    A%nrow = nrow
    A%ncol = ncol
    A%max_degree = 0

end subroutine cs_matrix_init



! Should we put a subroutine in here which sets the row/column oriented
! storage format for a matrix, or should we change the API so that the
! assembly routine now takes an additional argument "orientation" that
! determines the format?



!--------------------------------------------------------------------------!
subroutine cs_assemble(A,g)                                                !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    class(cs_graph), pointer, intent(in) :: g

    A%g => g

    ! Check whether the matrix is stored in row- or column-oriented format.
    ! If we take a CSR matrix A and use the CSC matvec operation
    !       A%csc_matvec(x,y),
    ! we fortuitously get an output vector y: y = transpose(A)*x.
    ! Likewise, CSR_matvec = CSC_transpose_matvec.
    ! We exploit this property by having
    !       csr_matvec_t => csc_matvec,     csc_matvec_t => csr_matvec.
    ! This saves us from writing excess code.
    if (A%orientation=='col') then
        A%ncol = g%n
        A%nrow = g%m

        A%find_entry => csc_find_entry

        A%matvec_impl => csc_matvec
        A%matvec_t_impl => csr_matvec

        A%left_permute_impl => cs_matrix_permute_vals
        A%right_permute_impl => cs_matrix_permute_ptrs
    else
        A%nrow = g%n
        A%ncol = g%m

        A%find_entry => csr_find_entry

        A%matvec_impl => csr_matvec
        A%matvec_t_impl => csc_matvec

        A%left_permute_impl => cs_matrix_permute_ptrs
        A%right_permute_impl => cs_matrix_permute_vals
    endif
    A%nnz = g%ne
    A%max_degree = g%max_degree

    allocate(A%val(A%nnz))
    A%val = 0.0_dp

end subroutine cs_assemble



!--------------------------------------------------------------------------!
subroutine cs_matrix_neighbors(A,i,nbrs)                                   !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(in)  :: i
    integer, intent(out) :: nbrs(:)

    nbrs = 0
    call A%g%neighbors(i,nbrs)

end subroutine cs_matrix_neighbors



!--------------------------------------------------------------------------!
function cs_get_value(A,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: cs_get_value
    ! local variables
    integer :: k

    cs_get_value = 0_dp
    k = A%find_entry(i,j)
    if (k/=-1) cs_get_value = A%val(k)

end function cs_get_value



!--------------------------------------------------------------------------!
subroutine cs_set_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%find_entry(i,j)
    if (k/=-1) then
        A%val(k) = val
    else
        call A%cs_set_value_not_preallocated(i,j,val)
    endif

end subroutine cs_set_value



!--------------------------------------------------------------------------!
subroutine cs_add_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%find_entry(i,j)
    if (k/=-1) then
        A%val(k) = A%val(k)+val
    else
        call A%cs_set_value_not_preallocated(i,j,val)
    endif

end subroutine cs_add_value



!--------------------------------------------------------------------------!
subroutine cs_sub_matrix_add(A,B)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    class(cs_matrix), intent(in)    :: B
    ! local variables
    integer :: i,j,k,indx,nbrs(B%max_degree)

    ! Jesus, really need to fix this thing
    do i=1,B%g%n
        do k=B%g%ptr(i),B%g%ptr(i+1)-1
            j = B%g%node(k)
            indx = A%g%find_edge(i,j)
            A%val(indx) = A%val(indx)+B%val(k)
        enddo
    enddo

end subroutine cs_sub_matrix_add



!--------------------------------------------------------------------------!
subroutine cs_left_permute(A,p)                                            !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_permute_impl(p)

end subroutine cs_left_permute



!--------------------------------------------------------------------------!
subroutine cs_right_permute(A,p)                                           !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_permute_impl(p)

end subroutine cs_right_permute



!--------------------------------------------------------------------------!
subroutine cs_matvec(A,x,y)                                                !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_impl(x,y)

end subroutine cs_matvec



!--------------------------------------------------------------------------!
subroutine cs_matvec_t(A,x,y)                                              !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_t_impl(x,y)

end subroutine cs_matvec_t



!--------------------------------------------------------------------------!
subroutine cs_set_value_not_preallocated(A,i,j,val)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    real(dp) :: val_temp(A%nnz)
    integer :: k

    if (A%orientation=='col') then
        call A%g%add_edge(j,i)
    else
        call A%g%add_edge(i,j)
    endif

    k = A%find_entry(i,j)

    val_temp = A%val
    deallocate(A%val)
    allocate(A%val(A%nnz+1))
    A%val(1:k-1) = val_temp(1:k-1)
    A%val(k) = val
    A%val(k+1:A%nnz+1) = val_temp(k:A%nnz)
    A%nnz = A%nnz+1
    A%max_degree = A%g%max_degree

end subroutine cs_set_value_not_preallocated





!==========================================================================!
!==========================================================================!
!== Implementations of matrix operations                                 ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
function csr_find_entry(A,i,j)                                             !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    integer :: csr_find_entry

    csr_find_entry = A%g%find_edge(i,j)

end function csr_find_entry



!--------------------------------------------------------------------------!
function csc_find_entry(A,i,j)                                             !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    integer :: csc_find_entry

    csc_find_entry = A%g%find_edge(j,i)

end function csc_find_entry



!--------------------------------------------------------------------------!
subroutine csr_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z

    do i=1,A%g%n
        z = 0_dp
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
            z = z+A%val(k)*x(j)
        enddo
        y(i) = z
    enddo

end subroutine csr_matvec



!--------------------------------------------------------------------------!
subroutine csc_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z

    y = 0_dp
    do i=1,A%g%n
        z = x(i)
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
            y(j) = y(j)+A%val(k)*z
        enddo
    enddo

end subroutine csc_matvec



!--------------------------------------------------------------------------!
subroutine cs_matrix_permute_vals(A,p)                                     !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine cs_matrix_permute_vals



!--------------------------------------------------------------------------!
subroutine cs_matrix_permute_ptrs(A,p)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,k,ptr(A%g%n+1)
    real(dp) :: val(A%nnz)

    do i=1,A%g%n
        ptr(p(i)+1) = A%g%ptr(i+1)-A%g%ptr(i)
    enddo

    ptr(1) = 1
    do i=1,A%g%n
        ptr(i+1) = ptr(i+1)+ptr(i)
    enddo

    do i=1,A%g%n
        do k=0,A%g%ptr(i+1)-A%g%ptr(i)-1
            val( ptr(p(i))+k ) = A%val( A%g%ptr(i)+k )
        enddo
    enddo

    A%val = val

    call A%g%left_permute(p)

end subroutine cs_matrix_permute_ptrs





end module cs_matrices

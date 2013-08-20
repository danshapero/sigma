module bcsr_matrices

use sparse_matrices
use block_sparse_matrices
use cs_graphs

implicit none


!--------------------------------------------------------------------------!
type, extends(block_sparse_matrix) :: bcsr_matrix                          !
!--------------------------------------------------------------------------!
    class(cs_graph), pointer :: g
    real(dp), allocatable :: val(:,:,:)
contains
    procedure :: init => bcsr_init
    procedure :: assemble => bcsr_assemble
    procedure :: neighbors => bcsr_matrix_neighbors
    procedure :: get_value => bcsr_get_value
    procedure :: set_value => bcsr_set_value, add_value => bcsr_add_value
    procedure :: get_block => bcsr_get_block
    procedure :: set_block => bcsr_set_block, add_block => bcsr_add_block
    procedure :: sub_matrix_add => bcsr_sub_matrix_add
    procedure :: left_permute => bcsr_left_permute, &
                & right_permute => bcsr_right_permute
    procedure :: matvec => bcsr_matvec, matvec_t => bcsr_matvec_t
    procedure :: block_matvec => bcsr_block_matvec, &
                & block_matvec_t => bcsr_block_matvec_t
    procedure :: l_block_matvec => bcsr_l_block_matvec, &
                & l_block_matvec_t => bcsr_l_block_matvec_t
    procedure :: r_block_matvec => bcsr_r_block_matvec, &
                & r_block_matvec_t => bcsr_r_block_matvec_t
    !procedure, private :: bcsr_set_block_not_preallocated
end type bcsr_matrix


contains



!--------------------------------------------------------------------------!
subroutine bcsr_init(A,nrow,ncol)                                          !
!--------------------------------------------------------------------------!
    class(bcsr_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol

    A%nrow = nrow
    A%ncol = ncol
    A%max_degree = 0

end subroutine bcsr_init



!--------------------------------------------------------------------------!
subroutine bcsr_assemble(A,g)                                              !
!--------------------------------------------------------------------------!
    class(bcsr_matrix), intent(inout) :: A
    class(cs_graph), pointer, intent(in) :: g

    A%g => g

    A%nr = A%nrow/g%n
    A%nc = A%ncol/g%m
    A%nnz = g%ne
    A%max_degree = g%max_degree

    allocate(A%val(A%nr,A%nc,A%nnz))
    A%val = 0.0_dp

end subroutine bcsr_assemble



!--------------------------------------------------------------------------!
subroutine bcsr_matrix_neighbors(A,i,nbrs)                                 !
!--------------------------------------------------------------------------!
    class(bcsr_matrix), intent(in) :: A
    integer, intent(in)  :: i
    integer, intent(out) :: nbrs(:)

    nbrs = 0
    call A%g%neighbors(i,nbrs)

end subroutine bcsr_matrix_neighbors



!--------------------------------------------------------------------------!
function bcsr_get_value(A,i,j)                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: bcsr_get_value
    ! local variables
    integer :: k,rb,cb,ib,jb

    rb = (i-1)/A%nr+1
    cb = (j-1)/A%nc+1
    ib = i-A%nr*rb
    jb = j-A%nc*cb

    bcsr_get_value = 0_dp
    k = A%g%find_edge(rb,cb)
    if (k/=-1) bcsr_get_value = A%val(ib,jb,k)

end function bcsr_get_value



!--------------------------------------------------------------------------!
subroutine bcsr_set_value(A,i,j,val)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,rb,cb,ib,jb

    rb = (i-1)/A%nr+1
    cb = (j-1)/A%nc+1
    ib = i-A%nr*rb
    jb = j-A%nc*cb

    k = A%g%find_edge(rb,cb)
    if (k/=-1) then
        A%val(ib,jb,k) = val
    endif
    ! Change this to an else statement and add soemthing to add values in
    ! space which hasn't been preallocated yet

end subroutine bcsr_set_value



!--------------------------------------------------------------------------!
subroutine bcsr_add_value(A,i,j,val)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,rb,cb,ib,jb

    rb = (i-1)/A%nr+1
    cb = (j-1)/A%nc+1
    ib = i-A%nr*rb
    jb = j-A%nc*cb

    k = A%g%find_edge(rb,cb)
    if (k/=-1) then
        A%val(ib,jb,k) = A%val(ib,jb,k)+val
    endif
    ! Change this to an else statement and add soemthing to add values in
    ! space which hasn't been preallocated yet

end subroutine bcsr_add_value



!--------------------------------------------------------------------------!
subroutine bcsr_get_block(A,i,j,vals)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp), intent(out) :: vals(:,:)
    ! local variables
    integer :: k

    k = A%g%find_edge(i,j)
    vals = A%val(:,:,k)

end subroutine bcsr_get_block



!--------------------------------------------------------------------------!
subroutine bcsr_set_block(A,i,j,vals)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: vals(:,:)
    ! local variables
    integer :: k

    k = A%g%find_edge(i,j)
    A%val(:,:,k) = vals

end subroutine bcsr_set_block



!--------------------------------------------------------------------------!
subroutine bcsr_add_block(A,i,j,vals)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: vals(:,:)
    ! local variables
    integer :: k

    k = A%g%find_edge(i,j)
    A%val(:,:,k) = A%val(:,:,k)+vals

end subroutine bcsr_add_block



!--------------------------------------------------------------------------!
subroutine bcsr_sub_matrix_add(A,B)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(inout)   :: A
    class(bcsr_matrix), intent(in)      :: B
    ! local variables
    integer :: i,j,k,indx

    do i=1,B%nrow
        do k=B%g%ia(i),B%g%ia(i+1)-1
            j = B%g%ja(k)
            indx = A%g%find_edge(i,j)
            A%val(:,:,indx) = A%val(:,:,indx)+B%val(:,:,k)
        enddo
    enddo

end subroutine bcsr_sub_matrix_add



!--------------------------------------------------------------------------!
subroutine bcsr_left_permute(A,p)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,k,ia(A%g%n+1)
    real(dp) :: val(A%nnz,A%nr,A%nc)

    do i=1,A%g%n
        ia(p(i)+1) = A%g%ia(i+1)-A%g%ia(i)
    enddo

    ia(1) = 1
    do i=1,A%g%n
        ia(i+1) = ia(i+1)+ia(i)
    enddo

    do i=1,A%g%n
        do k=0,A%g%ia(i+1)-A%g%ia(i)-1
            val(ia(p(i))+k,:,:) = A%val(A%g%ia(i)+k,:,:)
        enddo
    enddo

    A%val = val

    call A%g%left_permute(p)

end subroutine bcsr_left_permute



!--------------------------------------------------------------------------!
subroutine bcsr_right_permute(A,p)                                         !
!--------------------------------------------------------------------------!
    class(bcsr_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine bcsr_right_permute



!--------------------------------------------------------------------------!
subroutine bcsr_matvec(A,x,y)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(A%nr)

    do i=1,A%g%n
        z = 0.0_dp
        do k=A%g%ia(i),A%g%ia(i+1)-1
            j = A%g%ja(k)
            do l=1,A%nc
                z = z+A%val(:,l,k)*x(A%nc*(j-1)+1)
            enddo
        enddo
        y(A%nr*(i-1)+1:A%nr*i) = z
    enddo

end subroutine bcsr_matvec



!--------------------------------------------------------------------------!
subroutine bcsr_matvec_t(A,x,y)                                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(A%nr)

    do i=1,A%g%n
        z = x(A%nr*(i-1)+1:A%nr*i)
        do k=A%g%ia(i),A%g%ia(i+1)-1
            j = A%g%ja(k)
            do l=1,A%nc
                y(A%nc*(j-1)+l) = y(A%nc*(j-1)+1) &
                    & +dot_product(A%val(:,k,l),z)
            enddo
        enddo
    enddo

end subroutine bcsr_matvec_t



!--------------------------------------------------------------------------!
subroutine bcsr_block_matvec(A,x,y)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:,:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(A%nr)

    do i=1,A%g%n
        z = 0.0_dp
        do k=A%g%ia(i),A%g%ia(i+1)-1
            j = A%g%ja(k)
            do l=1,A%nc
                z = z+A%val(:,l,k)*x(l,j)
            enddo
        enddo
        y(:,i) = z
    enddo

end subroutine bcsr_block_matvec



!--------------------------------------------------------------------------!
subroutine bcsr_block_matvec_t(A,x,y)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:,:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(A%nc)

    do i=1,A%g%n
        z = x(:,i)
        do k=A%g%ia(i),A%g%ia(i+1)-1
            j = A%g%ja(k)
            do l=1,A%nc
                y(l,j) = y(l,j)+dot_product(A%val(:,l,k),z)
            enddo
        enddo
    enddo

end subroutine bcsr_block_matvec_t



!--------------------------------------------------------------------------!
subroutine bcsr_l_block_matvec(A,x,y)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:,:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(A%nr)

    do i=1,A%g%n
        z = 0.0_dp
        do k=A%g%ia(i),A%g%ia(i+1)-1
            j = A%g%ja(k)
            do l=1,A%nc
                z = z+A%val(:,l,k)*x(A%nc*(j-1)+l)
            enddo
        enddo
        y(:,i) = z
    enddo

end subroutine bcsr_l_block_matvec



!--------------------------------------------------------------------------!
subroutine bcsr_l_block_matvec_t(A,x,y)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:,:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(A%nr)

    do i=1,A%g%n
        z = x(A%nr*(i-1)+1:A%nr*i)
        do k=A%g%ia(i),A%g%ia(i+1)-1
            j = A%g%ja(k)
            do l=1,A%nc
                y(l,j) = y(l,j)+dot_product(A%val(:,l,k),z)
            enddo
        enddo
    enddo

end subroutine bcsr_l_block_matvec_t



!--------------------------------------------------------------------------!
subroutine bcsr_r_block_matvec(A,x,y)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(A%nr)

    do i=1,A%g%n
        z = 0.0_dp
        do k=A%g%ia(i),A%g%ia(i+1)-1
            j = A%g%ja(k)
            do l=1,A%nc
                z = z+A%val(:,l,k)*x(l,j)
            enddo
        enddo
        y(A%nr*(i-1)+1:A%nr*i) = z
    enddo

end subroutine bcsr_r_block_matvec



!--------------------------------------------------------------------------!
subroutine bcsr_r_block_matvec_t(A,x,y)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcsr_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(A%nr)

    do i=1,A%g%n
        z = x(:,i)
        do k=A%g%ia(i),A%g%ia(i+1)-1
            j = A%g%ja(k)
            do l=1,A%nc
                y(A%nc*(j-1)+l) = y(A%nc*(j-1)+l) &
                    & +dot_product(A%val(:,l,k),z)
            enddo
        enddo
    enddo

end subroutine bcsr_r_block_matvec_t






end module bcsr_matrices

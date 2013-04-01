!--------------------------------------------------------------------------!
module linalg_mod                                                          !
!--------------------------------------------------------------------------!

    use abstract_matrix_mod
    use workspace_mod

    implicit none



type :: linalg_work_arrays
    integer :: nn
    real(kind(1d0)), allocatable :: p(:),q(:),r(:),z(:)
contains
    procedure :: init_work_arrays
end type linalg_work_arrays


contains



!--------------------------------------------------------------------------!
subroutine init_work_arrays(work,nn)                                       !
!--------------------------------------------------------------------------!
    implicit none
    integer, intent(in) :: nn
    class (linalg_work_arrays), intent(inout) :: work
    allocate( work%p(nn),work%q(nn),work%r(nn),work%z(nn) )
    work%p = 0.d0
    work%q = 0.d0
    work%r = 0.d0
    work%z = 0.d0

end subroutine init_work_arrays



!--------------------------------------------------------------------------!
subroutine cg(A,x,b,tol,work)                                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    real(kind(1d0)), intent(in) :: tol
    type(workspace), intent(inout) :: work
    ! local variables
    real(kind(1d0)) :: alpha,beta,res2,dpr
    integer :: ind
    real(kind(1d0)), pointer :: p(:),q(:),r(:),z(:)

    ind = work%find_free_work_vectors(4)
    p => work%w(:, ind )
    q => work%w(:,ind+1)
    r => work%w(:,ind+2)
    z => work%w(:,ind+3)

    call A%matvec(x,q)
    r = b-q
    p = r
    res2 = dot_product(r,r)

    do while ( dsqrt(res2)>tol )
        call A%matvec(p,q)
        dpr = dot_product(p,q)
        alpha = res2/dpr
        x   = x+alpha*p
        r = r-alpha*q
        dpr = dot_product(r,r)
        beta = dpr/res2
        p = r+beta*p
        res2 = dpr
    enddo

end subroutine cg



!--------------------------------------------------------------------------!
subroutine subblock_cg(A,x,b,tol,work,i1,i2)                               !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    real(kind(1d0)), intent(in) :: tol
    class(workspace), intent(inout) :: work
    integer, intent(in) :: i1,i2
    ! local variables
    real(kind(1d0)) :: alpha,beta,res2,dpr
    integer :: i,ind
    real(kind(1d0)), pointer :: p(:),q(:),r(:),z(:)

    ind = work%find_free_work_vectors(4)
    p => work%w(:, ind )
    q => work%w(:,ind+1)
    r => work%w(:,ind+2)
    z => work%w(:,ind+3)

    call A%subblock_matvec(x,q,i1,i2,i1,i2)
    r = b-q
    do i=1,i1-1
        r(i) = 0.d0
    enddo
    do i=i2+1,size(x)
        r(i) = 0.d0
    enddo
    p = r
    res2 = dot_product(r,r)

    do while ( dsqrt(res2)>tol )
        call A%subblock_matvec(p,q,i1,i2,i1,i2)
        dpr = dot_product(p,q)
        alpha = res2/dpr
        x = x+alpha*p
        r = r-alpha*q
        dpr = dot_product(r,r)
        beta = dpr/res2
        p = r+beta*p
        res2 = dpr
    enddo

end subroutine subblock_cg



!--------------------------------------------------------------------------!
subroutine subset_cg(A,x,b,tol,work,setlist,set)                           !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    real(kind(1d0)), intent(in) :: tol
    class(workspace), intent(inout) :: work
    integer, dimension(:), intent(in) :: setlist
    integer, intent(in) :: set
    ! local variables
    real(kind(1d0)) :: alpha,beta,res2,dpr
    integer :: i,ind
    real(kind(1d0)), pointer :: p(:),q(:),r(:),z(:)

    ind = work%find_free_work_vectors(4)
    p => work%w(:, ind )
    q => work%w(:,ind+1)
    r => work%w(:,ind+2)
    z => work%w(:,ind+3)

    call A%subset_matvec(x,q,setlist,set,set)
    r = b-q
    do i=1,size(x)
        if (setlist(i) /= set) r(i) = 0.d0
    enddo
    p = r
    res2 = dot_product(r,r)

    do while( dsqrt(res2)>tol )
        call A%subset_matvec(p,q,setlist,set,set)
        dpr = dot_product(p,q)
        alpha = res2/dpr
        x = x+alpha*p
        r = r-alpha*q
        dpr = dot_product(r,r)
        beta = dpr/res2
        p = r+beta*p
        res2 = dpr
    enddo

end subroutine subset_cg



end module linalg_mod

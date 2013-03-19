!--------------------------------------------------------------------------!
module linalg_mod                                                          !
!--------------------------------------------------------------------------!

    use abstract_matrix_mod

    implicit none
    real(kind(1d0)), allocatable, private :: p(:),q(:),r(:),z(:)


contains



!--------------------------------------------------------------------------!
subroutine allocate_linalg_workarrays(nn)                                  !
!--------------------------------------------------------------------------!
    implicit none
    integer, intent(in) :: nn

    allocate( p(nn),q(nn),r(nn),z(nn) )
    p = 0.d0
    q = 0.d0
    r = 0.d0
    z = 0.d0

end subroutine allocate_linalg_workarrays


!--------------------------------------------------------------------------!
subroutine cg(A,x,b,tol)                                                   !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (sparse_matrix), intent(in) :: A
    real(kind(1d0)), dimension(:), intent(inout) :: x
    real(kind(1d0)), dimension(:), intent(in) :: b
    real(kind(1d0)), intent(in) :: tol
    ! local variables
    real(kind(1d0)) :: alpha,beta,res2,dpr

    call A%matvec(x,q)
    r = b-q
    p = r
    res2 = dot_product(r,r)

    do while ( dsqrt(res2)>tol )
        call A%matvec(p,q)
        dpr = dot_product(p,q)
        alpha = res2/dpr
        x = x+alpha*p
        r = r-alpha*q
        dpr = dot_product(r,r)
        beta = dpr/res2
        p = r+beta*p
        res2 = dpr
    enddo

end subroutine cg



!--------------------------------------------------------------------------!
subroutine subblock_cg(A,x,b,tol,i1,i2)                                    !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (sparse_matrix), intent(in) :: A
    real(kind(1d0)), dimension(:), intent(inout) :: x
    real(kind(1d0)), dimension(:), intent(in) :: b
    real(kind(1d0)), intent(in) :: tol
    integer, intent(in) :: i1,i2
    ! local variables
    real(kind(1d0)) :: alpha,beta,res2,dpr
    integer :: i
    
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
subroutine subset_cg(A,x,b,tol,setlist,set)                                !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (sparse_matrix), intent(in) :: A
    real(kind(1d0)), dimension(:), intent(inout) :: x
    real(kind(1d0)), dimension(:), intent(in) :: b
    real(kind(1d0)), intent(in) :: tol
    integer, dimension(:), intent(in) :: setlist
    integer, intent(in) :: set
    ! local variables
    real(kind(1d0)) :: alpha,beta,res2,dpr
    integer :: i

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

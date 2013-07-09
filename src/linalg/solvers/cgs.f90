module cgs

    use matrix
    use solver

    implicit none




type, extends(iterative_solver) :: cgs_solver
    real(kind(1d0)), allocatable :: p(:),q(:),r(:),r0(:),z(:),u(:)
contains
    procedure :: init => cgs_init
    procedure :: solve => cgs_solve
    procedure :: destroy => cgs_destroy
end type cgs_solver



contains


!--------------------------------------------------------------------------!
subroutine cgs_init(solver,nn,tolerance)                                   !
!--------------------------------------------------------------------------!
    implicit none
    class(cgs_solver), intent(inout) :: solver
    integer, intent(in) :: nn
    real(kind(1d0)), intent(in) :: tolerance

    solver%nn = nn
    solver%tolerance = tolerance
    solver%iterations = 0

    allocate(solver%p(nn),solver%q(nn),solver%r(nn),solver%r0(nn), &
        & solver%z(nn),solver%u(nn))

end subroutine cgs_init



!--------------------------------------------------------------------------!
subroutine cgs_solve(solver,A,x,b,pc,mask)                                 !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(cgs_solver), intent(inout) :: solver
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    class(preconditioner), intent(inout) :: pc
    integer, intent(in) :: mask(:)
    ! local variables
    real(kind(1d0)) :: alpha,beta,dpr,res2

    associate( p=>solver%p, q=>solver%q, r=>solver%r, r0=>solver%r0, &
        z=>solver%z, u=>solver%u )

    call A%matvec(x,z)
    r = b-z
    r(mask) = 0.d0
    r0 = r
    p = r0
    u = r0
    res2 = dot_product(r,r)
    do while (dsqrt(res2)>solver%tolerance)
        dpr = dot_product(r,r0)
        call A%matvec(p,z)
        z(mask) = 0.d0
        alpha = dpr/dot_product(z,r0)
        q = u-alpha*z
        x = x+alpha*(u+q)
        call A%matvec(u+q,z)
        z(mask) = 0.d0
        r = r-alpha*z
        beta = dot_product(r,r0)/dpr
        u = r+beta*q
        p = u+beta*(q+beta*p)
        res2 = dot_product(r,r)
        solver%iterations = solver%iterations+1
    enddo

    end associate

end subroutine cgs_solve



!--------------------------------------------------------------------------!
subroutine cgs_destroy(solver)                                             !
!--------------------------------------------------------------------------!
    implicit none
    class(cgs_solver), intent(inout) :: solver

    deallocate( solver%p,solver%q,solver%r,solver%r0,solver%u,solver%z )
    solver%iterations = 0
    solver%nn = 0
    solver%tolerance = 0.d0

end subroutine cgs_destroy




end module cgs

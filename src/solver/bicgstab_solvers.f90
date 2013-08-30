module bicgstab_solvers

use types
use sparse_matrices
use iterative_solvers

implicit none


!--------------------------------------------------------------------------!
type, extends(iterative_solver) :: bicgstab_solver                         !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: p(:), q(:), r(:), r0(:), v(:), s(:), t(:), z(:)
    real(dp) :: alpha, beta, omega, rho, rho_old
contains
    procedure :: init => bicgstab_init
    procedure :: pc_solve => bicgstab_pc_solve
    procedure :: no_pc_solve => bicgstab_no_pc_solve
    procedure :: free => bicgstab_free
end type bicgstab_solver

contains



!--------------------------------------------------------------------------!
subroutine bicgstab_init(solver,nn,tolerance)                              !
!--------------------------------------------------------------------------!
    class(bicgstab_solver), intent(inout) :: solver
    integer, intent(in)                   :: nn
    real(dp), intent(in), optional        :: tolerance

    solver%nn = nn
    if (present(tolerance)) then
        solver%tolerance = tolerance
    else
        solver%tolerance = 1.0d-8
    endif

    allocate( solver%p(nn), solver%q(nn), solver%r(nn), solver%r0(nn), &
        & solver%v(nn), solver%s(nn), solver%t(nn), solver%z(nn) )

end subroutine bicgstab_init



!--------------------------------------------------------------------------!
subroutine bicgstab_pc_solve(solver,A,x,b,pc)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bicgstab_solver), intent(inout) :: solver
    class(sparse_matrix), intent(in)      :: A
    real(dp), intent(inout)               :: x(:)
    real(dp), intent(in)                  :: b(:)
    class(preconditioner), intent(inout)  :: pc
    ! local variables
    real(dp) :: res2

    associate( p => solver%p, q => solver%q, r => solver%r, &
        & r0 => solver%r0, v => solver%v, s => solver%s, t => solver%t, &
        & z => solver%z, alpha => solver%alpha, rho => solver%rho, &
        & rho_old => solver%rho_old, beta => solver%beta, &
        & omega => solver%omega )

    call A%matvec(x,q)
    z = b-q
    call pc%precondition(A,r0,z)
    r = r0

    rho = 1.0_dp
    rho_old = 1.0_dp
    alpha = 1.0_dp
    omega = 1.0_dp

    v = 0.0_dp
    p = 0.0_dp

    res2 = dot_product(r,r)

    do while ( dsqrt(res2)>solver%tolerance )
        rho = dot_product(r0,r)
        beta = rho/rho_old*alpha/omega
        p = r+beta*(p-omega*v)
        call A%matvec(p,z)
        call pc%precondition(A,v,z)

        alpha = rho/dot_product(r0,v)
        s = r-alpha*v
        call A%matvec(s,z)
        call pc%precondition(A,t,z)
        omega = dot_product(s,t)/dot_product(t,t)
        x = x+alpha*p+omega*s
        r = s-omega*t

        rho_old = rho
        res2 = dot_product(r,r)
    enddo

    end associate

end subroutine bicgstab_pc_solve



!--------------------------------------------------------------------------!
subroutine bicgstab_no_pc_solve(solver,A,x,b)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bicgstab_solver), intent(inout) :: solver
    class(sparse_matrix), intent(in)      :: A
    real(dp), intent(inout)               :: x(:)
    real(dp), intent(in)                  :: b(:)
    ! local variables
    real(dp) :: res2

    associate( p => solver%p, q => solver%q, r => solver%r, &
        & r0 => solver%r0, v => solver%v, s => solver%s, t => solver%t, &
        & z => solver%z, alpha => solver%alpha, rho => solver%rho, &
        & rho_old => solver%rho_old, beta => solver%beta, &
        & omega => solver%omega )

    call A%matvec(x,q)
    r0 = b-q
    r = r0

    rho = 1.0_dp
    rho_old = 1.0_dp
    alpha = 1.0_dp
    omega = 1.0_dp

    v = 0.0_dp
    p = 0.0_dp

    res2 = dot_product(r,r)

    do while ( dsqrt(res2)>solver%tolerance )
        rho = dot_product(r0,r)
        beta = rho/rho_old*alpha/omega
        p = r+beta*(p-omega*v)

        call A%matvec(p,v)
        alpha = rho/dot_product(r0,v)
        s = r-alpha*v

        call A%matvec(s,t)
        omega = dot_product(s,t)/dot_product(t,t)
        if (isnan(omega)) omega = 0.0_dp
        x = x+alpha*p+omega*s
        r = s-omega*t

        res2 = dot_product(r,r)
        rho_old = rho
    enddo

    end associate

end subroutine bicgstab_no_pc_solve



!--------------------------------------------------------------------------!
subroutine bicgstab_free(solver)                                           !
!--------------------------------------------------------------------------!
    class(bicgstab_solver), intent(inout) :: solver

    deallocate(solver%p, solver%q, solver%r, solver%r0, solver%s, &
        & solver%t, solver%v, solver %z)
    solver%alpha = 0.0_dp
    solver%rho = 0.0_dp
    solver%rho_old = 0.0_dp
    solver%omega = 0.0_dp

end subroutine bicgstab_free



end module bicgstab_solvers

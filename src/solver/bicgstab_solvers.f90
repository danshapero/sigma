module bicgstab_solvers

use types, only: dp
use linear_operator_interface

implicit none


!--------------------------------------------------------------------------!
type, extends(linear_solver) :: bicgstab_solver                            !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: p(:), q(:), r(:), r0(:), v(:), s(:), t(:), z(:)
    real(dp) :: alpha, beta, omega, rho, rho_old
contains
    procedure :: init => bicgstab_init
    procedure :: linear_solve => bicgstab_solve
    procedure :: linear_solve_pc => bicgstab_solve_pc
    procedure :: free => bicgstab_free
end type bicgstab_solver

contains



!--------------------------------------------------------------------------!
subroutine bicgstab_init(solver,A)                                         !
!--------------------------------------------------------------------------!
    class(bicgstab_solver), intent(inout) :: solver
    class(linear_operator), intent(in) :: A

    if (A%nrow/=A%ncol) then
        print *, 'Cannot make a BiCG solver for a non-square matrix.'
        print *, 'Terminating.'
        call exit(1)
    endif

    solver%nn = A%nrow
    solver%tolerance = 1.0d-8


    allocate( solver%p(solver%nn), solver%q(solver%nn), &
        & solver%r(solver%nn), solver%r0(solver%nn), &
        & solver%v(solver%nn), solver%s(solver%nn), &
        & solver%t(solver%nn), solver%z(solver%nn) )

end subroutine bicgstab_init



!--------------------------------------------------------------------------!
subroutine bicgstab_solve(solver,A,x,b)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bicgstab_solver), intent(inout) :: solver
    class(linear_operator), intent(in)    :: A
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

end subroutine bicgstab_solve



!--------------------------------------------------------------------------!
subroutine bicgstab_solve_pc(solver,A,x,b,pc)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bicgstab_solver), intent(inout) :: solver
    class(linear_operator), intent(in)    :: A
    real(dp), intent(inout)               :: x(:)
    real(dp), intent(in)                  :: b(:)
    class(linear_solver), intent(inout)   :: pc
    ! local variables
    real(dp) :: res2

    associate( p => solver%p, q => solver%q, r => solver%r, &
        & r0 => solver%r0, v => solver%v, s => solver%s, t => solver%t, &
        & z => solver%z, alpha => solver%alpha, rho => solver%rho, &
        & rho_old => solver%rho_old, beta => solver%beta, &
        & omega => solver%omega )

    call A%matvec(x,q)
    z = b-q
    call pc%solve(A,r0,z)
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
        call pc%solve(A,v,z)

        alpha = rho/dot_product(r0,v)
        s = r-alpha*v
        call A%matvec(s,z)
        call pc%solve(A,t,z)
        omega = dot_product(s,t)/dot_product(t,t)
        x = x+alpha*p+omega*s
        r = s-omega*t

        rho_old = rho
        res2 = dot_product(r,r)
    enddo

    end associate

end subroutine bicgstab_solve_pc



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

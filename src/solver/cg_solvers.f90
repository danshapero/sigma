module cg_solvers

use types, only: dp
use linear_operator_interface

implicit none


!--------------------------------------------------------------------------!
type, extends(linear_solver) :: cg_solver                                  !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: p(:), q(:), r(:), z(:)
    integer :: iterations
contains
    procedure :: init => cg_init
    procedure :: linear_solve => cg_solve
    procedure :: linear_solve_pc => cg_solve_pc
    procedure :: free => cg_free
end type cg_solver


contains



!--------------------------------------------------------------------------!
subroutine cg_init(solver,A)                                               !
!--------------------------------------------------------------------------!
    class(cg_solver), intent(inout) :: solver
    class(linear_operator), intent(in) :: A

    if (A%nrow/=A%ncol) then
        print *, 'Cannot make a CG solver for a non-square matrix.'
        print *, 'Terminating.'
        call exit(1)
    endif

    solver%nn = A%nrow
    solver%tolerance = 1.0d-8
    solver%iterations = 0

    allocate( solver%p(solver%nn), solver%q(solver%nn), &
        & solver%r(solver%nn), solver%z(solver%nn) )

    solver%p = 0.0_dp
    solver%q = 0.0_dp
    solver%r = 0.0_dp
    solver%z = 0.0_dp

end subroutine cg_init



!--------------------------------------------------------------------------!
subroutine cg_solve(solver,A,x,b)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cg_solver), intent(inout)    :: solver
    class(linear_operator), intent(in) :: A
    real(dp), intent(inout)            :: x(:)
    real(dp), intent(in)               :: b(:)
    ! local variables
    real(dp) :: alpha, beta, res2, dpr

    associate( p=>solver%p, q=>solver%q, r=>solver%r )

    call A%matvec(x,q)
    r = b-q
    p = r
    res2 = dot_product(r,r)

    do while( dsqrt(res2)>solver%tolerance )
        call A%matvec(p,q)
        dpr = dot_product(p,q)
        alpha = res2/dpr
        x = x+alpha*p
        r = r-alpha*q

        dpr = dot_product(r,r)
        beta = dpr/res2
        p = r+beta*p
        res2 = dpr

        solver%iterations = solver%iterations+1
    enddo

    end associate

end subroutine cg_solve



!--------------------------------------------------------------------------!
subroutine cg_solve_pc(solver,A,x,b,pc)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cg_solver), intent(inout)     :: solver
    class(linear_operator), intent(in)  :: A
    real(dp), intent(inout)             :: x(:)
    real(dp), intent(in)                :: b(:)
    class(linear_solver), intent(inout) :: pc
    ! local variables
    real(dp) :: alpha, beta, dpr, res2

    associate( p=>solver%p, q=>solver%q, r=>solver%r, z=>solver%z )

    z = x
    call A%matvec(z,q)
    r = b-q
    call pc%solve(A,z,r)
    p = z
    res2 = dot_product(r,z)

    do while( dsqrt(res2)>solver%tolerance )
        call A%matvec(p,q)
        dpr = dot_product(p,q)
        alpha = res2/dpr
        x = x+alpha*p
        r = r-alpha*q

        call pc%solve(A,z,r)

        dpr = dot_product(r,z)
        beta = dpr/res2
        p = z+beta*p
        res2 = dpr

        solver%iterations = solver%iterations+1
    enddo

    end associate

end subroutine cg_solve_pc



!--------------------------------------------------------------------------!
subroutine cg_free(solver)                                                 !
!--------------------------------------------------------------------------!
    class(cg_solver), intent(inout) :: solver

    solver%nn = 0
    solver%iterations = 0
    solver%tolerance = 0.0_dp

    deallocate( solver%p, solver%q, solver%r, solver%z )

end subroutine cg_free





end module cg_solvers

module cg_solvers

use types
use sparse_matrices
use iterative_solvers

implicit none


!--------------------------------------------------------------------------!
type, extends(iterative_solver) :: cg_solver                               !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: p(:), q(:), r(:), z(:)
contains
    procedure :: init => cg_init
    procedure :: pc_solve => cg_pc_solve
    procedure :: no_pc_solve => cg_no_pc_solve
    procedure :: free => cg_free
end type cg_solver


contains



!--------------------------------------------------------------------------!
subroutine cg_init(solver,nn,tolerance)                                    !
!--------------------------------------------------------------------------!
    class(cg_solver), intent(inout) :: solver
    integer, intent(in)             :: nn
    real(dp), intent(in), optional  :: tolerance

    solver%nn = nn
    if (present(tolerance)) then
        solver%tolerance = tolerance
    else
        solver%tolerance = 1.0d-8
    endif

    allocate( solver%p(nn), solver%q(nn), solver%r(nn), solver%z(nn) )

end subroutine cg_init



!--------------------------------------------------------------------------!
subroutine cg_pc_solve(solver,A,x,b,pc)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cg_solver), intent(inout)         :: solver
    class(sparse_matrix), intent(in)        :: A
    real(dp), intent(inout)                 :: x(:)
    real(dp), intent(in)                    :: b(:)
    class(preconditioner), intent(inout)    :: pc
    ! local variables
    real(dp) :: alpha, beta, dpr, res2

    associate( p=>solver%p, q=>solver%q, r=>solver%r, z=>solver%z )

    z = x
    call A%matvec(z,q)
    r = b-q
    call pc%precondition(A,z,r)
    p = z
    res2 = dot_product(r,z)

    do while( dsqrt(res2)>solver%tolerance )
        call A%matvec(p,q)
        dpr = dot_product(p,q)
        alpha = res2/dpr
        x = x+alpha*p
        r = r-alpha*q

        call pc%precondition(A,z,r)

        dpr = dot_product(r,z)
        beta = dpr/res2
        p = z+beta*p
        res2 = dpr

        solver%iterations = solver%iterations+1
    enddo

    end associate

end subroutine cg_pc_solve



!--------------------------------------------------------------------------!
subroutine cg_no_pc_solve(solver,A,x,b)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cg_solver), intent(inout)     :: solver
    class(sparse_matrix), intent(in)    :: A
    real(dp), intent(inout)             :: x(:)
    real(dp), intent(in)                :: b(:)
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

end subroutine cg_no_pc_solve



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



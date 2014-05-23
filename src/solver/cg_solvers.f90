module cg_solvers

use types, only: dp
use linear_operator_interface

implicit none


!--------------------------------------------------------------------------!
type, extends(linear_solver) :: cg_solver                                  !
!--------------------------------------------------------------------------!
    ! scratch vectors for CG iteration
    real(dp), allocatable :: p(:), q(:), r(:), z(:)
    integer :: iterations

    ! parameters determining solver behavior
    real(dp) :: tolerance
    logical, private :: params_set = .false.
contains
    ! Methods requiredby linear solver interface
    procedure :: setup => cg_setup
    procedure :: linear_solve => cg_solve
    procedure :: linear_solve_pc => cg_solve_pc
    procedure :: destroy => cg_destroy

    ! Methods specific to CG solvers
    procedure :: set_params => cg_set_params
end type cg_solver


contains



!--------------------------------------------------------------------------!
function cg(tolerance)                                                     !
!--------------------------------------------------------------------------!
    real(dp), intent(in), optional :: tolerance
    class(linear_solver), pointer :: cg

    allocate(cg_solver::cg)
    select type(cg)
        type is(cg_solver)
            call cg%set_params(tolerance)
    end select

end function cg



!--------------------------------------------------------------------------!
subroutine cg_setup(solver,A)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cg_solver), intent(inout) :: solver
    class(linear_operator), intent(in) :: A
    ! local variables
    integer :: nn

    ! Error handling
    if (A%ncol/=A%nrow) then
        print *, 'Cannot make a CG solver for a non-square matrix'
        print *, 'Terminating.'
        call exit(1)
    endif

    ! Set the matrix dimension for this solver
    nn = A%nrow
    solver%nn = nn

    ! Set the iteration count to 0
    solver%iterations = 0

    ! If the solver parameters have not been set, call the set_params
    ! subroutine to set them all to their defaults
    if (.not.solver%params_set) call solver%set_params()

    ! If the solver hasn't been initialized yet, allocate the work vectors
    if (.not.solver%initialized) then
        allocate( solver%p(nn), solver%q(nn), solver%r(nn), solver%z(nn) )
        solver%initialized = .true.
    endif

    ! Zero out all the work vectors
    solver%p = 0.0_dp
    solver%q = 0.0_dp
    solver%r = 0.0_dp
    solver%z = 0.0_dp

end subroutine cg_setup



!--------------------------------------------------------------------------!
subroutine cg_set_params(solver,tolerance)                                 !
!--------------------------------------------------------------------------!
    class(cg_solver), intent(inout) :: solver
    real(dp), intent(in), optional :: tolerance

    ! If the user has specified a tolerance for the iterative solver, then
    ! set the solver's tolerance accordingly,
    if (present(tolerance)) then
        solver%tolerance = tolerance
    else
    ! Otherwise, set the tolerance to 1.0e-16.
        solver%tolerance = 1.0d-16
    endif

    solver%params_set = .true.

end subroutine cg_set_params



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
subroutine cg_destroy(solver)                                              !
!--------------------------------------------------------------------------!
    class(cg_solver), intent(inout) :: solver

    solver%nn = 0
    solver%tolerance = 0.0_dp
    solver%iterations = 0

    deallocate( solver%p, solver%q, solver%r, solver%z )

    solver%initialized = .false.
    solver%params_set = .false.

end subroutine cg_destroy





end module cg_solvers

module gcr_solvers

use types, only: dp
use linear_operator_interface

implicit none


!--------------------------------------------------------------------------!
type, extends(linear_solver) :: gcr_solver                                 !
!--------------------------------------------------------------------------!
    ! scratch vectors for GCR iteration
    real(dp), allocatable :: ap(:), p(:), q(:), r(:), z(:)
    integer :: iterations

    ! parameters determining solver behavior
    real(dp) :: tolerance
    logical, private :: params_set = .false.
contains
    ! Methods required by linear solver itnerface
    procedure :: setup => gcr_setup
    procedure :: linear_solve => gcr_solve
    procedure :: linear_solve_pc => gcr_solve_pc
    procedure :: destroy => gcr_destroy

    ! Methods specific to GCR solvers
    procedure :: set_params => gcr_set_params
end type gcr_solver


contains



!--------------------------------------------------------------------------!
function gcr(tolerance)                                                    !
!--------------------------------------------------------------------------!
    real(dp), intent(in), optional :: tolerance
    class(linear_solver), pointer :: gcr

    allocate(gcr_solver::gcr)
    select type(gcr)
        type is(gcr_solver)
            call gcr%set_params(tolerance)
    end select

end function gcr



!--------------------------------------------------------------------------!
subroutine gcr_setup(solver,A)                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(gcr_solver), intent(inout) :: solver
    class(linear_operator), intent(in) :: A
    ! local variables
    integer :: nn

    ! Error handling
    if (A%ncol/=A%nrow) then
        print *, 'Cannot make a GCR solver for a non-square matrix'
        print *, 'Terminating.'
        call exit(1)
    endif

    ! Set the iteration count to 0
    solver%iterations = 0

    ! If the solver parameters have not been set, call the set_params
    ! subroutine to set them all to their defaults
    if (.not.solver%params_set) call solver%set_params()

    ! If the solver hasn't been initialized yet, allocate the work vectors
    if (.not.solver%initialized) then
        allocate( solver%ap(nn), solver%p(nn), solver%q(nn), &
            & solver%r(nn), solver%z(nn) )
        solver%initialized = .true.
    endif

    ! Zero out all the work vectors
    solver%ap = 0.0_dp
    solver%p = 0.0_dp
    solver%q = 0.0_dp
    solver%r = 0.0_dp
    solver%z = 0.0_dp

end subroutine gcr_setup



!--------------------------------------------------------------------------!
subroutine gcr_set_params(solver,tolerance)                                !
!--------------------------------------------------------------------------!
    class(gcr_solver), intent(inout) :: solver
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

end subroutine gcr_set_params



!--------------------------------------------------------------------------!
subroutine gcr_solve(solver,A,x,b)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(gcr_solver), intent(inout)    :: solver
    class(linear_operator), intent(in)  :: A
    real(dp), intent(inout)             :: x(:)
    real(dp), intent(in)                :: b(:)
    ! local variables
    real(dp) :: alpha, beta, res2, dpr

    associate( ap=>solver%ap, p=>solver%p, q=>solver%p, r=>solver%r)

    call A%matvec(x,q)
    r = b-q
    p = r
    res2 = dot_product(r,r)

    ! q will always store A*r
    call A%matvec(r,q)
    ap = q

    do while(dsqrt(res2)>solver%tolerance)
        dpr = dot_product(r,q)
        alpha = dpr/dot_product(ap,ap)
        x = x+alpha*p
        r = r-alpha*q

        call A%matvec(r,q)
        beta = dot_product(r,q)/dpr

        p = r+beta*p
        ap = q+beta*ap

        res2 = dot_product(r,r)

        solver%iterations = solver%iterations+1
    enddo

    end associate

end subroutine gcr_solve



!--------------------------------------------------------------------------!
subroutine gcr_solve_pc(solver,A,x,b,pc)                                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(gcr_solver), intent(inout)    :: solver
    class(linear_operator), intent(in)  :: A
    real(dp), intent(inout)             :: x(:)
    real(dp), intent(in)                :: b(:)
    class(linear_solver), intent(inout) :: pc
    ! local variables
    real(dp) :: alpha, beta, res2, dpr

    associate( ap=>solver%ap, p=>solver%p, q=>solver%p, r=>solver%r, &
        & z=>solver%z)

    call A%matvec(x,q)
    z = b-q
    call pc%solve(A,r,q)
    p = r

    ! q will always store A*r
    call A%matvec(r,q)
    ap = q

    ! z will always store M^{-1}*A*p
    call pc%solve(A,z,ap)

    do while(dsqrt(res2)>solver%tolerance)
        dpr = dot_product(r,q)
        alpha = dpr/dot_product(ap,z)
        x = x+alpha*p
        r = r-alpha*z

        call A%matvec(r,q)
        beta = dot_product(r,q)/dpr

        p = r+beta*p
        ap = q+beta*ap

        call pc%solve(A,z,ap)

        res2 = dot_product(r,r)

        solver%iterations = solver%iterations+1
    enddo

    end associate

end subroutine gcr_solve_pc



!--------------------------------------------------------------------------!
subroutine gcr_destroy(solver)                                             !
!--------------------------------------------------------------------------!
    class(gcr_solver), intent(inout) :: solver

    solver%nn = 0
    solver%tolerance = 0.0_dp
    solver%iterations = 0

    deallocate( solver%ap, solver%p, solver%q, solver%r, solver%z )

    solver%initialized = .false.
    solver%params_set = .false.

end subroutine gcr_destroy




end module gcr_solvers

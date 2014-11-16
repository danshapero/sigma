module bicgstab_solvers

use types, only: dp
use linear_operator_interface

implicit none


!--------------------------------------------------------------------------!
type, extends(linear_solver) :: bicgstab_solver                            !
!--------------------------------------------------------------------------!
    ! scratch vectors and scalars for BiCG algorithm
    real(dp), allocatable :: p(:), q(:), r(:), r0(:), v(:), s(:), t(:), z(:)
    real(dp) :: alpha, beta, omega, rho, rho_old
    integer :: iterations

    ! parameters determining the behavior of the solver
    real(dp) :: tolerance
    logical, private :: params_set = .false.
contains
    ! Methods required by the linear solver interface
    procedure :: setup => bicgstab_setup
    procedure :: linear_solve => bicgstab_solve
    procedure :: linear_solve_pc => bicgstab_solve_pc
    procedure :: destroy => bicgstab_destroy

    ! Methods specific to BiCG-Stab solvers
    procedure :: set_params => bicgstab_set_params
end type bicgstab_solver

contains



!--------------------------------------------------------------------------!
function bicgstab(tolerance)                                               !
!--------------------------------------------------------------------------!
    real(dp), intent(in), optional :: tolerance
    class(linear_solver), pointer :: bicgstab

    allocate(bicgstab_solver::bicgstab)
    select type(bicgstab)
        type is(bicgstab_solver)
            call bicgstab%set_params(tolerance)
    end select

end function bicgstab



!--------------------------------------------------------------------------!
subroutine bicgstab_setup(solver,A)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bicgstab_solver), intent(inout) :: solver
    clasS(linear_operator), intent(in) :: A
    ! local variables
    integer :: nn

    ! Error handling
    if (A%ncol/=A%nrow) then
        print *, 'Cannot make a BiCG-Stab solver for non-square matrix.'
        print *, 'Terminating.'
        call exit(1)
    endif

    ! Set the matrix dimension for this solver
    nn = A%nrow
    solver%nn = nn

    ! Set the iteration count to 0
    solver%iterations = 0

    ! If the solver parameters have not been set, call the set_params
    ! subroutine to set them all to defaults
    if (.not.solver%params_set) call solver%set_params()

    ! If the solver hasn't yet been initialized, allocate the work vectors
    if (.not.solver%initialized) then
        allocate( solver%p(nn), solver%q(nn), &
            & solver%r(nn), solver%r0(nn), &
            & solver%v(nn), solver%s(nn), &
            & solver%t(nn), solver%z(nn) )

        solver%initialized = .true.
    endif

    ! Zero out the work vectors
    solver%p = 0.0_dp
    solver%q = 0.0_dp
    solver%r = 0.0_dp
    solver%r0 = 0.0_dp
    solver%v = 0.0_dp
    solver%s = 0.0_dp
    solver%t = 0.0_dp
    solver%z = 0.0_dp

end subroutine bicgstab_setup



!--------------------------------------------------------------------------!
subroutine bicgstab_set_params(solver,tolerance)                           !
!--------------------------------------------------------------------------!
    class(bicgstab_solver), intent(inout) :: solver
    real(dp), intent(in), optional :: tolerance

    ! If the user has specified a tolerane for the iterative solver, then
    ! set the solver's tolerane accordingly,
    if (present(tolerance)) then
        solver%tolerance = tolerance
    else
    ! Otherwise, set the tolerance to 1.0e-16.
        solver%tolerance = 1.0d-16
    endif

    solver%params_set = .true.

end subroutine bicgstab_set_params



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

        solver%iterations = solver%iterations + 1
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

        solver%iterations = solver%iterations + 1
    enddo

    end associate

end subroutine bicgstab_solve_pc



!--------------------------------------------------------------------------!
subroutine bicgstab_destroy(solver)                                        !
!--------------------------------------------------------------------------!
    class(bicgstab_solver), intent(inout) :: solver

    deallocate(solver%p, solver%q, solver%r, solver%r0, solver%s, &
        & solver%t, solver%v, solver %z)
    solver%alpha = 0.0_dp
    solver%rho = 0.0_dp
    solver%rho_old = 0.0_dp
    solver%omega = 0.0_dp

    solver%initialized = .false.
    solver%params_set = .false.

end subroutine bicgstab_destroy



end module bicgstab_solvers

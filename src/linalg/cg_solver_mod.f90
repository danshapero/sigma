module cg_solver_mod

    use sparse_matrix_mod
    use iterative_solver_mod

    implicit none




type, extends(iterative_solver) :: cg_solver
    real(kind(1d0)), allocatable :: p(:),q(:),r(:),z(:)
contains
    procedure :: init => cg_init
    procedure :: solve => cg_solve
    procedure :: destroy => cg_destroy
end type cg_solver



contains


!--------------------------------------------------------------------------!
subroutine cg_init(solver,nn,tolerance)                                    !
!--------------------------------------------------------------------------!
    implicit none
    class(cg_solver), intent(inout) :: solver
    integer, intent(in) :: nn
    real(kind(1d0)), intent(in) :: tolerance

    solver%nn = nn
    solver%tolerance = tolerance
    solver%iterations = 0

    allocate(solver%p(nn),solver%q(nn),solver%r(nn),solver%z(nn))

end subroutine cg_init



!--------------------------------------------------------------------------!
subroutine cg_solve(solver,A,x,b,pc,mask)                                  !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(cg_solver), intent(inout) :: solver
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    class(preconditioner), intent(inout) :: pc
    integer, intent(in) :: mask(:)
    ! local variables
    real(kind(1d0)) :: alpha,beta,dpr,res2

    associate( p=>solver%p, q=>solver%q, r=>solver%r, z=>solver%z )

    z = x
    z(mask) = 0.d0
    call A%matvec(z,q)
    r = b-q
    r(mask) = 0.d0
    call pc%precondition(A,z,r,mask)
    p = z
    res2 = dot_product(r,z)

    do while ( dsqrt(res2)>solver%tolerance )
        call A%matvec(p,q)
        q(mask) = 0.d0
        dpr = dot_product(p,q)
        alpha = res2/dpr
        x = x+alpha*p
        r = r-alpha*q
        call pc%precondition(A,z,r,mask)
        dpr = dot_product(r,z)
        beta = dpr/res2
        p = z+beta*p
        res2 = dpr
        solver%iterations = solver%iterations+1
    enddo

    end associate

end subroutine cg_solve



!--------------------------------------------------------------------------!
subroutine cg_destroy(solver)                                              !
!--------------------------------------------------------------------------!
    implicit none
    class(cg_solver), intent(inout) :: solver

    deallocate( solver%p, solver%q, solver%r, solver%z )
    solver%iterations = 0
    solver%nn = 0
    solver%tolerance = 0.d0

end subroutine cg_destroy




end module cg_solver_mod

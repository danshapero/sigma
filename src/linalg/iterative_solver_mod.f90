module iterative_solver_mod

    use sparse_matrix_mod

    implicit none
    integer, save :: null_mask(0)


!--------------------------------------------------------------------------!
! Iterative solver abstract data type                                      !
!--------------------------------------------------------------------------!
type, abstract :: iterative_solver
    integer :: nn,iterations
    real(kind(1d0)) :: tolerance
contains
    procedure(solver_init_interface), deferred :: init
    procedure(solve_interface), deferred :: solve
    procedure(solver_destroy_interface), deferred :: destroy
end type iterative_solver



!--------------------------------------------------------------------------!
! Preconditioner abstract data type                                        !
!--------------------------------------------------------------------------!
type, abstract :: preconditioner
    integer :: nn,level
contains
    procedure(preconditioner_init_interface), deferred :: init
    procedure(precondition_interface), deferred :: precondition
    procedure(preconditioner_clear_interface), deferred :: clear
end type preconditioner


! Null preconditioner
type, extends(preconditioner) :: nopc

contains
    procedure :: init => nopc_init
    procedure :: precondition => nopc_precondition
    procedure :: clear => nopc_clear
end type nopc



!--------------------------------------------------------------------------!
! Interfaces for all of solver and preconditioner type-bound procedures    !
!--------------------------------------------------------------------------!
abstract interface

!-------------------
! Solver interfaces
subroutine solver_init_interface(solver,nn,tolerance)
    import :: iterative_solver
    class(iterative_solver), intent(inout) :: solver
    integer, intent(in) :: nn
    real(kind(1d0)), intent(in) :: tolerance
end subroutine solver_init_interface

subroutine solve_interface(solver,A,x,b,pc,mask)
    import :: iterative_solver, sparse_matrix, preconditioner
    class(iterative_solver), intent(inout) :: solver
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    class(preconditioner), intent(inout) :: pc
    integer, intent(in) :: mask(:)
end subroutine solve_interface

subroutine solver_destroy_interface(solver)
    import :: iterative_solver
    class(iterative_solver), intent(inout) :: solver
end subroutine


!---------------------------
! Preconditioner interfaces
subroutine preconditioner_init_interface(pc,A,level)
    import :: preconditioner, sparse_matrix
    class(preconditioner), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: level
end subroutine preconditioner_init_interface

subroutine precondition_interface(pc,A,x,b,mask)
    import :: preconditioner, sparse_matrix
    class(preconditioner), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: mask(:)
end subroutine precondition_interface

subroutine preconditioner_clear_interface(pc)
    import :: preconditioner
    class(preconditioner), intent(inout) :: pc
end subroutine preconditioner_clear_interface


end interface



contains



!--------------------------------------------------------------------------!
subroutine nopc_init(pc,A,level)                                           !
!--------------------------------------------------------------------------!
    implicit none
    class(nopc), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: level

    pc%nn = A%nrow
    pc%level = level

end subroutine nopc_init



!--------------------------------------------------------------------------!
subroutine nopc_precondition(pc,A,x,b,mask)                                !
!--------------------------------------------------------------------------!
    implicit none
    class(nopc), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: mask(:)

    x = b

end subroutine nopc_precondition



!--------------------------------------------------------------------------!
subroutine nopc_clear(pc)                                                  !
!--------------------------------------------------------------------------!
    implicit none
    class(nopc), intent(inout) :: pc

    ! Do nothing

end subroutine nopc_clear



end module iterative_solver_mod

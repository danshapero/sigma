module iterative_solvers

use types
use sparse_matrices

implicit none



!--------------------------------------------------------------------------!
type, abstract :: iterative_solver                                         !
!--------------------------------------------------------------------------!
    integer :: nn, iterations
    real(dp) :: tolerance
contains
    procedure(init_solver_ifc), deferred    :: init
    procedure(pc_solve_ifc), deferred       :: pc_solve
    procedure(no_pc_solve_ifc), deferred    :: no_pc_solve
    procedure(free_solver_ifc), deferred    :: free
    generic :: solve => pc_solve, no_pc_solve
end type iterative_solver



!--------------------------------------------------------------------------!
type, abstract :: preconditioner                                           !
!--------------------------------------------------------------------------!
    integer :: nn, level
contains
    procedure(init_preconditioner_ifc), deferred    :: init
    procedure(precondition_ifc), deferred           :: precondition
    procedure(free_preconditioner_ifc), deferred    :: free
end type preconditioner


!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine init_solver_ifc(solver,nn,tolerance)
        import :: iterative_solver, dp
        class(iterative_solver), intent(inout)  :: solver
        integer, intent(in)                     :: nn
        real(dp), intent(in), optional          :: tolerance
    end subroutine init_solver_ifc

    subroutine pc_solve_ifc(solver,A,x,b,pc)
        import :: iterative_solver, sparse_matrix, preconditioner, dp
        class(iterative_solver), intent(inout)  :: solver
        class(sparse_matrix), intent(in)        :: A
        real(dp), intent(inout)                 :: x(:)
        real(dp), intent(in)                    :: b(:)
        class(preconditioner), intent(inout)    :: pc
    end subroutine pc_solve_ifc

    subroutine no_pc_solve_ifc(solver,A,x,b)
        import :: iterative_solver, sparse_matrix, dp
        class(iterative_solver), intent(inout)  :: solver
        class(sparse_matrix), intent(in)        :: A
        real(dp), intent(inout)                 :: x(:)
        real(dp), intent(in)                    :: b(:)
    end subroutine no_pc_solve_ifc

    subroutine free_solver_ifc(solver)
        import :: iterative_solver
        class(iterative_solver), intent(inout)  :: solver
    end subroutine free_solver_ifc


    subroutine init_preconditioner_ifc(pc,A,level)
        import :: preconditioner, sparse_matrix
        class(preconditioner), intent(inout)    :: pc
        class(sparse_matrix), intent(in)        :: A
        integer, intent(in), optional           :: level
    end subroutine init_preconditioner_ifc

    subroutine precondition_ifc(pc,A,x,b)
        import :: preconditioner, sparse_matrix, dp
        class(preconditioner), intent(inout)    :: pc
        class(sparse_matrix), intent(in)        :: A
        real(dp), intent(inout)                 :: x(:)
        real(dp), intent(in)                    :: b(:)
    end subroutine precondition_ifc

    subroutine free_preconditioner_ifc(pc)
        import :: preconditioner
        class(preconditioner), intent(inout)    :: pc
    end subroutine free_preconditioner_ifc

end interface




end module iterative_solvers

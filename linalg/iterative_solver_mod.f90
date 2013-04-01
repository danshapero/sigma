module iterative_solver_mod

    use sparse_matrix_mod

    implicit none



type, abstract :: iterative_solver
    integer :: nn,iterations
    real(kind(1d0)) :: tolerance
contains
    procedure(init_interface), deferred :: init
    procedure(solve_interface), deferred :: solve
    procedure(subblock_solve_interface), deferred :: subblock_solve
    procedure(subset_solve_interface), deferred :: subset_solve
    procedure(destroy_interface), deferred :: destroy
end type iterative_solver



abstract interface

subroutine init_interface(solver,nn,tolerance)
    import :: iterative_solver
    class(iterative_solver), intent(inout) :: solver
    integer, intent(in) :: nn
    real(kind(1d0)), intent(in) :: tolerance
end subroutine init_interface

subroutine solve_interface(solver,A,x,b)
    import :: iterative_solver, sparse_matrix
    class(iterative_solver), intent(inout) :: solver
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
end subroutine solve_interface

subroutine subblock_solve_interface(solver,A,x,b,i1,i2)
    import :: iterative_solver, sparse_matrix
    class(iterative_solver), intent(inout) :: solver
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: i1,i2
end subroutine subblock_solve_interface

subroutine subset_solve_interface(solver,A,x,b,setlist,set)
    import :: iterative_solver, sparse_matrix
    class(iterative_solver), intent(inout) :: solver
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: setlist(:),set
end subroutine subset_solve_interface

subroutine destroy_interface(solver)
    import :: iterative_solver
    class(iterative_solver), intent(inout) :: solver
end subroutine

end interface



contains
!-----------------
! Nothing! This module exists solely to define the iterative_solver type

end module iterative_solver_mod

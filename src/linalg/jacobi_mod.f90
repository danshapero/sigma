module jacobi_mod

    use sparse_matrix_mod
    use iterative_solver_mod

    implicit none




type, extends(preconditioner) :: jacobi
    real(kind(1d0)), allocatable, private :: diag(:)
contains
    procedure :: init => jacobi_init
    procedure :: precondition => jacobi_precondition
    procedure :: clear => jacobi_clear
end type jacobi



contains


!--------------------------------------------------------------------------!
subroutine jacobi_init(pc,A,level)                                         !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(jacobi), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: level
    ! local variables
    integer :: i

    pc%nn = A%nrow
    pc%level = level
    if (.not.allocated(pc%diag)) then
        allocate( pc%diag(pc%nn) )
    endif

    do i=1,A%nrow
        pc%diag(i) = A%get_value(i,i)
    enddo

end subroutine jacobi_init



!--------------------------------------------------------------------------!
subroutine jacobi_precondition(pc,A,x,b,mask)                              !
!--------------------------------------------------------------------------!
    implicit none
    class(jacobi), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: mask(:)

    x = b/pc%diag
    x(mask) = 0.d0

end subroutine jacobi_precondition



!--------------------------------------------------------------------------!
subroutine jacobi_clear(pc)                                                !
!--------------------------------------------------------------------------!
    implicit none
    class(jacobi), intent(inout) :: pc

    pc%diag = 0.d0
    pc%level = 0

end subroutine jacobi_clear



end module jacobi_mod

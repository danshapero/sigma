module nullpc_mod

    use sparse_matrix_mod
    use iterative_solver_mod

    implicit none




type, extends(preconditioner) :: nullpc

contains
    procedure :: init => nullpc_init
    procedure :: precondition => nullpc_precondition
    procedure :: subblock_precondition => nullpc_subblock_precondition
    procedure :: subset_precondition => nullpc_subset_precondition
    procedure :: clear => nullpc_clear
end type nullpc



contains



!--------------------------------------------------------------------------!
subroutine nullpc_init(pc,A,level)                                         !
!--------------------------------------------------------------------------!
    implicit none
    class(nullpc), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: level

    pc%nn = A%nrow
    pc%level = level

end subroutine nullpc_init



!--------------------------------------------------------------------------!
subroutine nullpc_precondition(pc,A,x,b)                                   !
!--------------------------------------------------------------------------!
    implicit none
    class(nullpc), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)

    x = b

end subroutine nullpc_precondition



!--------------------------------------------------------------------------!
subroutine nullpc_subblock_precondition(pc,A,x,b,i1,i2)                    !
!--------------------------------------------------------------------------!
    implicit none
    class(nullpc), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: i1,i2

    x = b

end subroutine nullpc_subblock_precondition



!--------------------------------------------------------------------------!
subroutine nullpc_subset_precondition(pc,A,x,b,setlist,set)                !
!--------------------------------------------------------------------------!
    implicit none
    class(nullpc), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: setlist(:),set

    x = b

end subroutine nullpc_subset_precondition



!--------------------------------------------------------------------------!
subroutine nullpc_clear(pc)                                                !
!--------------------------------------------------------------------------!
    implicit none
    class(nullpc), intent(inout) :: pc

    ! Do nothing. obviously.

end subroutine nullpc_clear



end module nullpc_mod

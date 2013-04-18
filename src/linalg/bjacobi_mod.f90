! This module does not work. We need a smarter system for picking out
! subblocks of the matrix, rather than just taking n x n chunks out in
! sequence. For example, do a breadth-first search and pick out n x n
! clusters which are then more likely to be connected.










module bjacobi_mod

    use sparse_matrix_mod
    use iterative_solver_mod

    implicit none




type, extends(preconditioner) :: bjacobi
    real(kind(1d0)), allocatable, private :: diag(:,:,:)
!    integer, allocatable, private :: ipiv(:,:)
    integer, private :: num_blocks
contains
    procedure :: init => bjacobi_init
    procedure :: precondition => bjacobi_precondition
    procedure :: clear => bjacobi_clear
end type bjacobi




contains


!--------------------------------------------------------------------------!
subroutine bjacobi_init(pc,A,level)                                        !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bjacobi), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: level
    ! local variables
    integer :: i,j,n,info,start,finish,rows(level)

    pc%nn = A%nrow
    pc%level = level
    pc%num_blocks = ceiling( pc%nn/real(level) )

    if (.not.allocated(pc%diag)) then
        allocate( pc%diag(level,level,pc%num_blocks) )
    endif
    pc%diag = 0.d0
 
    do i=1,pc%num_blocks
        start = level*(i-1)+1
        finish = min(level*i,pc%nn)
        n = finish-start+1
        do j=1,n
            rows(j) = start+j-1
        enddo
        pc%diag(1:n,1:n,i) = A%get_values( rows(1:n), rows(1:n) )

        call dpotrf('L',n,pc%diag(1:n,1:n,i),n,info)
        if (info/=0) print *, i,info

    enddo

end subroutine bjacobi_init



!--------------------------------------------------------------------------!
subroutine bjacobi_precondition(pc,A,x,b,mask)                             !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bjacobi), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: mask(:)
    ! local variables
    integer :: i,n,start,finish,info

    x = b
    x(mask) = 0.d0

    do i=1,pc%num_blocks
        start = pc%level*(i-1)+1
        finish = min(pc%level*i,pc%nn)
        n = finish-start+1
        call dpotrs('L',n,1,pc%diag(1:n,1:n,i),n,x(start:finish),n,info)
    enddo

    x(mask) = 0.d0

end subroutine bjacobi_precondition



!--------------------------------------------------------------------------!
subroutine bjacobi_clear(pc)                                               !
!--------------------------------------------------------------------------!
    implicit none
    class(bjacobi), intent(inout) :: pc

    pc%level = 0
    pc%diag = 0.d0

end subroutine bjacobi_clear






end module bjacobi_mod
    

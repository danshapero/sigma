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
    procedure :: subblock_precondition => bjacobi_subblock_precondition
    procedure :: subset_precondition => bjacobi_subset_precondition

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

    allocate( pc%diag(level,level,pc%num_blocks) )
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
subroutine bjacobi_precondition(pc,A,x,b)                                  !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bjacobi), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    ! local variables
    integer :: i,n,start,finish,info

    x = b

    do i=1,pc%num_blocks
        start = pc%level*(i-1)+1
        finish = min(pc%level*i,pc%nn)
        n = finish-start+1
        call dpotrs('L',n,1,pc%diag(1:n,1:n,i),n,x(start:finish),n,info)
    enddo

end subroutine bjacobi_precondition



!--------------------------------------------------------------------------!
subroutine bjacobi_subblock_precondition(pc,A,x,b,i1,i2)                   !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bjacobi), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: i1,i2
    ! local variables
    integer :: i,n,start,finish,si,fi,info
    real(kind(1d0)) :: v(pc%level)

    x = b

    do i=1,pc%num_blocks
        v = 0.d0

        start = pc%level*(i-1)+1
        finish = min(pc%level*i,pc%nn)
        n = finish-start+1
        si = max(start,i1)-start+1
        fi = min(finish,i2)-finish+n
        start = max(start,i1)
        finish = min(finish,i2)

        v(si:fi) = x(start:finish)

        if (n>0) then
            call dpotrs('L',n,1,pc%diag(1:n,1:n,i),n,v,n,info)
            x(start:finish) = v(si:fi)
        endif
    enddo

end subroutine bjacobi_subblock_precondition



!--------------------------------------------------------------------------!
subroutine bjacobi_subset_precondition(pc,A,x,b,setlist,set)               !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bjacobi), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: setlist(:), set
    ! local variables
    integer :: i,j,n,start,finish,info
    real(kind(1d0)) :: v(pc%level)

    x = b

    do i=1,pc%num_blocks
        v = 0.d0

        start = pc%level*(i-1)+1
        finish = min(pc%level*i,pc%nn)
        n = finish-start+1

        do j=0,n-1
            if (setlist(start+j) == set) v(1+j) = x(start+j)
        enddo

        if (n>0) then
            call dpotrs('L',n,1,pc%diag(1:n,1:n,i),n,v,n,info)
            do j=0,n-1
                if (setlist(start+j) == set) x(start+j) = v(1+j)
            enddo            
        endif
    enddo

end subroutine bjacobi_subset_precondition




end module bjacobi_mod
    

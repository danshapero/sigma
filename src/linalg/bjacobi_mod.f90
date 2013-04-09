module bjacobi_mod

    use sparse_matrix
    use iterative_solver_mod

    implicit none




type, extends(preconditioner) :: bjacobi
    real(kind(1d0)), allocatable, private :: diag(:,:,:)
    integer, allocatable, private :: ipiv(:,:)
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
    integer :: i,n,start,finish

    pc%nn = A%nrow
    pc%level = level
    pc%num_blocks = ceiling( pc%nn/real(level) )
    allocate( pc%diag(level,level,pc%num_blocks) )
    if (A%pos_def) then
        dns_solve => chol_dns_solve
    else
        allocate( ipiv(level,pc%num_blocks) )
        dns_solve => lu_dns_solve
    endif

    do i=1,pc%num_blocks
        start = level*(i-1)+1
        finish = min(level*i,pc%nn)
        n = finish-start+1
        pc%diag(1:n,1:n,i) = A%get_values( start:finish, start:finish )

        if (A%pos_def) then
            call dpotrf('L',n,pc%diag(1:n,1:n,i),n,n)
        else
            call dgetrf(n,n,pc%diag(1:n,1:n,i),n,ipiv(1:n,i),n)
        endif
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

    if (A%pos_def) then
        do i=1,pc%num_blocks
            start = pc%level*(i-1)+1
            finish = min(pc%level*i,pc%nn)
            n = finish-start+1
            call dpotrs('L',n,1,pc%diag(1:n,1:n,i),n,x(start:finish),n,info)
        enddo
    else
        do i=1,pc%num_blocks
            start = pc%level*(i-1)+1
            finish = min(pc%level*i,pc%nn)
            n = finish-start+1
            call dgetrs('N',n,1,pc%diag(1:n,1:n,i),n,pc%ipiv(1:n,i), &
                & x(start:finish),n,info)
        enddo
    endif

end subroutine bjacobi_precondition




end module bjacobi_mod
    

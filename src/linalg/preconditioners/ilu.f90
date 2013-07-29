module ilu

    use csr
    use matrix
    use solver

type, extends(preconditioner) :: ilu_preconditioner
    type(csr_matrix) :: LU
    real(kind(1d0)), allocatable :: D(:),q(:)
    logical :: pos_def
    logical :: initialized=.false.
contains
    procedure :: init => ilu_init
    procedure :: precondition => ilu_precondition
    procedure :: clear => ilu_clear
end type ilu_preconditioner



contains


!--------------------------------------------------------------------------!
subroutine ilu_init(pc,A,level)                                            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ilu_preconditioner), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: level
    ! local variables
    integer :: i,j,k,ptr1,ptr2,row,col
    integer, allocatable :: rows(:),cols(:)
    real(kind(1d0)) :: U
    real(kind(1d0)), allocatable :: vals(:)

    pc%nn = A%nrow
    pc%level = level

    pc%pos_def = A%pos_def

    allocate(rows(A%nnz),cols(A%nnz),vals(A%nnz))
    call A%convert_to_coo(rows,cols,vals)

    if (.not.pc%initialized) then
        allocate(pc%D(pc%nn),pc%q(pc%nn))
        call pc%LU%init(A%nrow,A%ncol,A%nnz,rows,cols)

        pc%initialized = .true.
    endif

    associate( LU=>pc%LU, D=>pc%D )

    do i=1,A%nnz
        call LU%set_value(rows(i),cols(i),vals(i))
    enddo

    deallocate(rows,cols,vals)

    do i=1,A%nrow
        do ptr1=LU%ia(i),LU%ia(i+1)-1
            k = LU%ja(ptr1)
            if ( k<i ) then
                LU%val(ptr1) = LU%val(ptr1)/D(k)
                do ptr2=ptr1+1,LU%ia(i+1)-1
                    j = LU%ja(ptr2)
                    U = LU%get_value(k,j)
                    LU%val(ptr2) = LU%val(ptr2)-LU%val(ptr1)*D(k)*U
                enddo
            elseif ( k==i ) then
                D(k) = LU%val(ptr1)
                LU%val(ptr1) = 1.d0
            else
                LU%val(ptr1) = LU%val(ptr1)/D(i)
            endif
        enddo
    enddo

    end associate

end subroutine ilu_init



!--------------------------------------------------------------------------!
subroutine ilu_precondition(pc,A,x,b,mask)                                 !
!--------------------------------------------------------------------------!
    implicit none
    class(ilu_preconditioner), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)
    integer, intent(in) :: mask(:)

    associate( LU => pc%LU, D => pc%D )

    x = b
    x(mask) = 0.d0
    call LU%forwardsolve(x)
    x(mask) = 0.d0
    x = x/D
    call LU%backsolve(x)
    x(mask) = 0.d0

    end associate

end subroutine ilu_precondition



!--------------------------------------------------------------------------!
subroutine ilu_clear(pc)                                                   !
!--------------------------------------------------------------------------!
    implicit none
    class(ilu_preconditioner), intent(inout) :: pc

    pc%D = 0.0
    pc%q = 0.d0
    pc%LU%val = 0.d0
    pc%pos_def = .false.

end subroutine ilu_clear


end module ilu

module ilu_mod

    use sparse_matrix_mod
    use iterative_solver_mod
    use csr_matrix_mod

type, extends(preconditioner) :: ilu
    class(sparse_matrix), allocatable :: LU
    real(kind(1d0)), allocatable :: D(:),q(:)
    logical :: pos_def
contains
    procedure :: init => ilu_init
    procedure :: precondition => ilu_precondition
    procedure :: subblock_precondition => ilu_subblock_precondition
    procedure :: subset_precondition => ilu_subset_precondition
    procedure :: clear => ilu_clear
end type ilu



contains


!--------------------------------------------------------------------------!
subroutine ilu_init(pc,A,level)                                            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(incomplete_cholesky), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: level
    ! local variables
    integer :: i,j,k,ptr1,ptr2,row,col,nbrs(A%max_degree)
    integer, allocatable :: rows(:),cols(:)
    real(kind(1d0)) :: U

    pc%nn = A%nrow
    pc%level = level

    pc%pos_def = A%pos_def

    allocate(csr_matrix::LU)
    allocate(D(pc%nn),q(pc%nn))

    type select(A)
        type is (csr_matrix)
            LU%nrow = A%nrow
            LU%ncol = A%ncol
            LU%nnz = A%nnz
            allocate(LU%ia(L%nrow+1),LU%ja(LU%nnz),LU%val(L%nnz))
            LU%ia = A%ia
            LU%ja = A%ja
            LU%val = A%val
        class default
            allocate(rows(A%nnz),cols(A%nnz),A%val)
            call A%convert_to_coo(rows,cols)
            LU%init(A%nrow,A%ncol,A%nnz,rows,cols,A%val)
            deallocate(rows,cols)
    end select

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
            elseif ( i==k ) then
                D(i) = LU%val(ptr1)
                LU%val(ptr1) = 1.d0
            else
                LU%val(ptr1) = LU%val(ptr1)/D(k)
            endif
        enddo
    enddo

end subroutine ilu_init



!--------------------------------------------------------------------------!
subroutine ilu_precondition(pc,A,x,b)                                      !
!--------------------------------------------------------------------------!
    implicit none
    class(ilu), intent(inout) :: pc
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    real(kind(1d0)), intent(in) :: b(:)

    call LU%forwardsolve(x,b)
    q = x/D
    call LU%backsolve(x,q)

end subroutine ilu_precondition



!--------------------------------------------------------------------------!
subroutine ilu_clear(pc)                                                   !
!--------------------------------------------------------------------------!
    implicit none
    class(ilu), intent(inout) :: pc

    deallocate(D,q)
    deallocate(LU)
    pc%pos_def = .false.

end subroutine ilu_clear


end module ilu_mod

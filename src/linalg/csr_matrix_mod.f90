module csr_matrix_mod

    use sparse_matrix_mod

    implicit none



type, extends(sparse_matrix) :: csr_matrix
    integer, allocatable :: ia(:), ja(:)
    real(kind(1d0)), allocatable :: val(:)
contains
    procedure :: init_matrix => csr_init_matrix
    procedure :: get_value => csr_get_value
    procedure :: get_values => csr_get_values
    procedure :: set_value => csr_set_value, add_value => csr_add_value
    procedure :: set_values => csr_set_values, add_values => csr_add_values
    procedure :: matvec => csr_matvec
    procedure :: subblock_matvec => csr_subblock_matvec
    procedure :: subset_matvec => csr_subset_matvec
    procedure :: subset_matrix_add => csr_subset_matrix_add
    procedure :: write_to_file => csr_write_to_file
end type csr_matrix



contains


!--------------------------------------------------------------------------!
subroutine csr_init_matrix(A,nrow,ncol,nnz,rows,cols)                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(inout) :: A
    integer, intent(in) :: nrow,ncol,nnz
    integer, intent(in), optional :: rows(:),cols(:)
    ! local variables
    integer :: i,j,k,n,startptr
    integer :: work(nrow)

    A%nrow = nrow
    A%ncol = ncol
    A%nnz = nnz

    allocate( A%ia(nrow+1), A%ja(nnz), A%val(nnz) )
    A%ia = 0
    A%ia(nrow+1) = nnz+1
    A%ja = 0
    A%val = 0.d0

    if (present(rows).and.present(cols)) then
        ! Ascertain how many non-zero entries are in each row
        work = 0
        do i=1,nnz
            work( rows(i) ) = work( rows(i) )+1
        enddo

        ! The maximum degree of the matrix is the greatest number of
        ! entries in any given row. (Strictly speaking not true for a 
        ! structurally asymmetric matrix, but for a CSR matrix this
        ! definition will do.)
        A%max_degree = maxval( work )

        ! Fill in the array ia; ia(i+1)-ia(i) = #of non-zeros in row i
        A%ia(1) = 1
        do i=1,nrow
            A%ia(i+1) = A%ia(i)+work(i)
        enddo

        ! Fill in the array ja; ja(j) = the column of non-zero entry #j
        work = 0
        do i=1,nnz
            startptr = A%ia( rows(i) )
            A%ja( startptr+work(rows(i)) ) = cols(i)
            work( rows(i) ) = work( rows(i) )+1
        enddo

        ! Sort the array ja so that ja(ia(i)),ja(ia(i)+1)-1 are in order
        do i=1,nrow
            do j=A%ia(i)+1,A%ia(i+1)-1
                k = j-1
                n = A%ja(j)
                do while( k>=A%ia(i) .and. A%ja(k)>n )
                    A%ja(k+1) = A%ja(k)
                    k = k-1
                enddo
                A%ja(k+1) = n
            enddo
        enddo

    endif

end subroutine csr_init_matrix






!==========================================================================!
!==========================================================================!
!===== Accessors & mutators                                            ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
real(kind(1d0)) function csr_get_value(A,i,j) result(value)                !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    ! local variables
    integer :: k

    do k=A%ia(i),A%ia(i+1)-1
        if ( A%ja(k) == j ) value = A%val(k)
    enddo

end function csr_get_value



!--------------------------------------------------------------------------!
function csr_get_values(A,rows,cols)                                       !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(in) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)) :: csr_get_values(size(rows),size(cols))
    ! local variables
    integer :: i,j,k

    do j=1,size(cols)
        do i=1,size(rows)
            do k=A%ia(rows(i)),A%ia(rows(i)+1)-1
                if ( A%ja(k)==cols(j) ) csr_get_values(i,j) = A%val(k)
            enddo
        enddo
    enddo
end function csr_get_values



!--------------------------------------------------------------------------!
subroutine csr_set_value(A,i,j,value)                                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: value
    ! local variables
    integer :: k

    do k=A%ia(i),A%ia(i+1)-1
        if ( A%ja(k) == j ) A%val(k) = value
    enddo

end subroutine csr_set_value



!--------------------------------------------------------------------------!
subroutine csr_add_value(A,i,j,value)                                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: value
    ! local variables
    integer :: k

    do k=A%ia(i),A%ia(i+1)-1
        if ( A%ja(k) == j ) A%val(k) = A%val(k)+value
    enddo

end subroutine csr_add_value



!--------------------------------------------------------------------------!
subroutine csr_set_values(A,rows,cols,values)                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: values(size(rows),size(cols))
    ! local variables
    integer :: i,j,k

    do j=1,size(cols)
        do i=1,size(rows)
            do k=A%ia(rows(i)),A%ia(rows(i)+1)-1
                if ( A%ja(k)==cols(j) ) A%val(k) = values(i,j)
            enddo
        enddo
    enddo

end subroutine csr_set_values



!--------------------------------------------------------------------------!
subroutine csr_add_values(A,rows,cols,values)                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: values(size(rows),size(cols))
    ! local variables
    integer :: i,j,k

    do j=1,size(cols)
        do i=1,size(rows)
            do k=A%ia(rows(i)),A%ia(rows(i)+1)-1
                if ( A%ja(k)==cols(j) ) A%val(k) = A%val(k)+values(i,j)
            enddo
        enddo
    enddo

end subroutine csr_add_values









!==========================================================================!
!==========================================================================!
!===== Matrix-vector multiplication routines                           ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine csr_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(in) :: A
    real(kind=8), intent(in) :: x(:)
    real(kind=8), intent(out) :: y(:)
    ! local variables
    real(kind=8) :: z
    integer :: i,j

    do i=1,A%nrow
        z = 0.d0
        do j=A%ia(i),A%ia(i+1)-1
            z = z+A%val(j)*x(A%ja(j))
        enddo
        y(i) = z
    enddo

end subroutine csr_matvec



!--------------------------------------------------------------------------!
subroutine csr_subblock_matvec(A,x,y,i1,i2,j1,j2)                          !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(in) :: A
    real(kind=8), intent(in) :: x(:)
    real(kind=8), intent(inout) :: y(:)
    integer, intent(in) :: i1,i2,j1,j2
    ! local variables
    real(kind=8) :: z
    integer :: i,j,k

    do i=i1,i2
        z = 0.d0
        do j=A%ia(i),A%ia(i+1)-1
            k = A%ja(j)
            if ( j1<=k .and. k<=j2 ) then
                z = z+A%val(j)*x(k)
            endif
        enddo
        y(i) = z
    enddo

end subroutine csr_subblock_matvec



!--------------------------------------------------------------------------!
subroutine csr_subset_matvec(A,x,y,setlist,set1,set2)                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(in) :: A
    real(kind=8), intent(in) :: x(:)
    real(kind=8), intent(inout) :: y(:)
    integer, dimension(:), intent(in) :: setlist
    integer, intent(in) :: set1,set2
    ! local variables
    real(kind=8) :: z
    integer :: i,j,k

    do i=1,A%nrow
        if ( setlist(i)==set1 ) then
            z = 0.d0
            do j=A%ia(i),A%ia(i+1)-1
                k = A%ja(j)
                if ( setlist(k)==set2 ) then
                    z = z+A%val(j)*x(k)
                endif
            enddo
            y(i) = z
        endif
    enddo

end subroutine csr_subset_matvec















!==========================================================================!
!==========================================================================!
!====== Matrix addition routines                                      =====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine csr_subset_matrix_add(A,B)                                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(inout) :: A
    class (csr_matrix), intent(in) :: B
    ! local variables
    integer :: i,j,k

    do i=1,A%nrow
        do j=A%ia(i),A%ia(i+1)-1
            do k=B%ia(i),B%ia(i+1)
                if ( B%ja(k)==A%ja(j) ) A%val(j) = A%val(j)+B%val(k)
            enddo
        enddo
    enddo

end subroutine csr_subset_matrix_add



!--------------------------------------------------------------------------!
subroutine csr_write_to_file(A,filename)                                   !
!--------------------------------------------------------------------------!
! Write out the data for a sparse matrix to a file in order to guarantee   !
! correctness                                                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(in) :: A
    character(len=*), intent(in) :: filename
    ! local variables
    integer :: n

    open(unit=10,file=trim(filename)//"ia.txt")
    do n=1,A%nrow+1
        write(10,*) A%ia(n)
    enddo
    close(10)

    open(unit=20,file=trim(filename)//"ja.txt")
    open(unit=30,file=trim(filename)//"a.txt")
    do n=1,A%nnz
        write(20,*) A%ja(n)
        write(30,*) A%val(n)
    enddo
    close(20)
    close(30)

end subroutine csr_write_to_file



end module csr_matrix_mod

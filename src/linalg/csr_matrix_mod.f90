module csr_matrix_mod

    use sparse_matrix_mod
    use omp_lib

    implicit none



type, extends(sparse_matrix) :: csr_matrix
    integer, allocatable :: ia(:), ja(:)
    real(kind(1d0)), allocatable :: val(:)
contains
    ! Constructor and accessors/mutators
    procedure :: init => csr_init
    procedure :: build => csr_build
    procedure :: get_value => csr_get_value
    procedure :: get_values => csr_get_values
    procedure :: get_neighbors => csr_get_neighbors
    procedure :: set_value => csr_set_value, add_value => csr_add_value
    procedure :: set_values => csr_set_values, add_values => csr_add_values
    procedure :: permute => csr_permute
    procedure :: subset_matrix_add => csr_subset_matrix_add
    ! matrix multiplication routines
    procedure :: matvec => csr_matvec
    ! forward- and back-solves for triangular systems
    procedure :: backsolve => csr_backsolve
    procedure :: forwardsolve => csr_forwardsolve
    ! routines for i/o and validation
    procedure :: convert_to_coo => csr_convert_to_coo
    procedure :: write_to_file => csr_write_to_file
    ! auxiliary routines
    procedure :: sort_ja
end type csr_matrix



contains







!==========================================================================!
!==========================================================================!
!===== Accessors & mutators                                            ====!
!==========================================================================!
!==========================================================================!



!--------------------------------------------------------------------------!
subroutine csr_init(A,nrow,ncol,nnz,rows,cols,vals,params)                 !
!--------------------------------------------------------------------------!
    implicit none
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: nrow,ncol,nnz
    integer, intent(in), optional :: rows(:),cols(:),params(:)
    real(kind(1d0)), intent(in), optional :: vals(:)

    A%nrow = nrow
    A%ncol = ncol
    A%nnz = nnz

    allocate( A%ia(nrow+1), A%ja(nnz), A%val(nnz) )
    A%ia = 0
    A%ia(nrow+1) = nnz+1
    A%ja = 0
    A%val = 0.d0

    if (present(rows).and.present(cols)) then
        call A%build(rows,cols,vals)
    endif

end subroutine csr_init



!--------------------------------------------------------------------------!
subroutine csr_build(A,rows,cols,vals)                                     !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in), optional :: vals(:)
    ! local variables
    integer :: i,j,k,n,startptr,nrow,ncol,nnz
    integer :: work(A%nrow)
    real(kind(1d0)) :: Aij

    nrow = A%nrow
    ncol = A%ncol
    nnz = A%nnz

    ! Fill in a work array; work(i) = # of non-zero entries in row i
    work = 0
    do i=1,nnz
        work( rows(i) ) = work( rows(i) )+1
    enddo

    A%max_degree = maxval(work)

    ! Fill in the array ia; ia(i+1)-ia(i) = # of non-zeros in row i
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

    if (present(vals)) then
        A%val = vals
    else
        A%val = 0.d0
    endif

    ! Sort the array ja: ja(ia(i)),ja(ia(i+1)-1) are in order
    call A%sort_ja()

end subroutine csr_build



!--------------------------------------------------------------------------!
function csr_get_value(A,i,j)                                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)) :: csr_get_value
    ! local variables
    integer :: k

    csr_get_value = 0.d0

    do k=A%ia(i),A%ia(i+1)-1
        if ( A%ja(k) == j ) csr_get_value = A%val(k)
    enddo

end function csr_get_value



!--------------------------------------------------------------------------!
function csr_get_values(A,rows,cols)                                       !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)) :: csr_get_values(size(rows),size(cols))
    ! local variables
    integer :: i,j,k

    csr_get_values = 0.d0

    do j=1,size(cols)
        do i=1,size(rows)
            do k=A%ia(rows(i)),A%ia(rows(i)+1)-1
                if ( A%ja(k)==cols(j) ) csr_get_values(i,j) = A%val(k)
            enddo
        enddo
    enddo

end function csr_get_values



!--------------------------------------------------------------------------!
function csr_get_neighbors(A,row)                                          !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    integer, intent(in) :: row
    integer :: csr_get_neighbors( A%max_degree )
    ! local variables
    integer :: start,finish

    csr_get_neighbors = 0
    start = A%ia(row)
    finish = A%ia(row+1)-1
    csr_get_neighbors(1:finish-start+1 ) = A%ja( start:finish )

end function csr_get_neighbors



!--------------------------------------------------------------------------!
subroutine csr_set_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: val
    ! local variables
    integer :: k

    do k=A%ia(i),A%ia(i+1)-1
        if ( A%ja(k) == j ) A%val(k) = val
    enddo

end subroutine csr_set_value



!--------------------------------------------------------------------------!
subroutine csr_add_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class (csr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: val
    ! local variables
    integer :: k

    do k=A%ia(i),A%ia(i+1)-1
        if ( A%ja(k) == j ) A%val(k) = A%val(k)+val
    enddo

end subroutine csr_add_value



!--------------------------------------------------------------------------!
subroutine csr_set_values(A,rows,cols,vals)                                !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: vals(size(rows),size(cols))
    ! local variables
    integer :: i,j,k

    do j=1,size(cols)
        do i=1,size(rows)
            do k=A%ia(rows(i)),A%ia(rows(i)+1)-1
                if ( A%ja(k)==cols(j) ) A%val(k) = vals(i,j)
            enddo
        enddo
    enddo

end subroutine csr_set_values



!--------------------------------------------------------------------------!
subroutine csr_add_values(A,rows,cols,vals)                                !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: vals(size(rows),size(cols))
    ! local variables
    integer :: i,j,k

    do j=1,size(cols)
        do i=1,size(rows)
            do k=A%ia(rows(i)),A%ia(rows(i)+1)-1
                if ( A%ja(k)==cols(j) ) A%val(k) = A%val(k)+vals(i,j)
            enddo
        enddo
    enddo

end subroutine csr_add_values



!--------------------------------------------------------------------------!
subroutine csr_permute(A,p)                                                !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    type(csr_matrix) :: B
    integer :: i,j,k,n,astart,bstart, q( A%nrow )
    real(kind(1d0)) :: Aij

    do i=1,A%nrow
        q( p(i) ) = i
    enddo

    allocate( B%ia( A%nrow+1 ), B%ja( A%nnz ), B%val( A%nnz ) )

    B%ia = A%ia
    B%ja = A%ja
    B%val = A%val

    A%ia = 0
    A%ja = 0
    A%val = 0.d0

    A%ia(1) = 1
    do i=1,A%nrow
        astart = A%ia(i)
        bstart = B%ia(q(i))
        A%ia(i+1) = astart+B%ia(q(i)+1)-B%ia(q(i))
        do j=0,A%ia(i+1)-A%ia(i)-1
            A%ja(astart+j) = p( B%ja(bstart+j) )
            A%val(astart+j) = B%val(bstart+j)
        enddo
    enddo

    ! Sort the array ja so that ja(ia(i)),ja(ia(i+1)-1) are in order
    call A%sort_ja()

end subroutine csr_permute



!--------------------------------------------------------------------------!
subroutine csr_subset_matrix_add(A,B)                                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    class(csr_matrix), intent(in) :: B
    ! local variables
    integer :: i,j,k

    !$omp parallel do private(j,k)
    do i=1,A%nrow
        do j=A%ia(i),A%ia(i+1)-1
            do k=B%ia(i),B%ia(i+1)-1
                if ( B%ja(k)==A%ja(j) ) A%val(j) = A%val(j)+B%val(k)
            enddo
        enddo
    enddo
    !$omp end parallel do

    A%symmetric = A%symmetric .and. B%symmetric
    A%pos_def = A%pos_def .and. B%pos_def
    A%diag_dominant = A%diag_dominant .and. B%diag_dominant

end subroutine csr_subset_matrix_add











!==========================================================================!
!==========================================================================!
!===== Matrix-vector multiplication routines                           ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine csr_matvec(A,x,y,rows,cols)                                     !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    real(kind(1d0)), intent(in) :: x(:)
    real(kind(1d0)), intent(out) :: y(:)
    integer, intent(in), optional :: rows(2),cols(2)
    ! local variables
    real(kind(1d0)) :: z
    integer :: i,j,ind,r(2),c(2)

    r = [1, A%nrow]
    c = [1, A%ncol]
    if (present(rows)) r = rows
    if (present(cols)) c = cols

    !$omp parallel do private(ind,j,z)
    do i=r(1),r(2)
        z = 0.d0
        do ind=A%ia(i),A%ia(i+1)-1
            j = A%ja(ind)
            if ( c(1)<=j .and. j<=c(2) ) z = z+A%val(ind)*x(j)
        enddo
        y(i) = z
    enddo
    !$omp end parallel do

end subroutine csr_matvec



!--------------------------------------------------------------------------!
subroutine csr_backsolve(A,x)                                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    ! local variables
    integer :: i,j
    real(kind(1d0)) :: Aii,z

    do i=A%nrow,1,-1
        z = x(i)
        do j=A%ia(i),A%ia(i+1)-1
            if (A%ja(j)>i) then
                z = z-A%val(j)*x(A%ja(j))
            elseif (A%ja(j)==i) then
                Aii = A%val(j)
            endif
        enddo
        x(i) = z/Aii
    enddo

end subroutine csr_backsolve



!--------------------------------------------------------------------------!
subroutine csr_forwardsolve(A,x)                                           !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    ! local variables
    integer :: i,j
    real(kind(1d0)) :: Aii,z

    do i=1,A%nrow
        z = x(i)
        do j=A%ia(i),A%ia(i+1)-1
            if (A%ja(j)<i) then
                z = z-A%val(j)*x(A%ja(j))
            elseif (A%ja(j)==i) then
                Aii = A%val(j)
            endif
        enddo
        x(i) = z/Aii
    enddo

end subroutine csr_forwardsolve









!==========================================================================!
!==========================================================================!
!======= i/o and validation routines                                   ====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine csr_convert_to_coo(A,rows,cols,vals)                            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    integer, intent(out) :: rows(:),cols(:)
    real(kind(1d0)), intent(out), optional :: vals(:)
    ! local variables
    integer :: i,j

    do i=1,A%nrow
        do j=A%ia(i),A%ia(i+1)-1
            rows(j) = i
        enddo
    enddo

    cols = A%ja
    if (present(vals)) vals = A%val

end subroutine csr_convert_to_coo



!--------------------------------------------------------------------------!
subroutine csr_write_to_file(A,filename)                                   !
!--------------------------------------------------------------------------!
! Write out the data for a sparse matrix to a file in order to guarantee   !
! correctness                                                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(in) :: A
    character(len=*), intent(in) :: filename
    ! local variables
    integer :: i,j

    open(unit=100,file=trim(filename)//".txt")
    write(100,1) A%nrow,A%nnz
1   format(' ',I5,'   ',I5)
    do i=1,A%nrow
        write(100,10) i
10      format(' ',I5,'  ')
        do j=A%ia(i),A%ia(i+1)-1
            write(100,20) A%ja(j),A%val(j)
20          format('     ',I5,'     ',g12.6)
        enddo
    enddo
    close(100)


end subroutine csr_write_to_file





!==========================================================================!
!==========================================================================!
!==== Auxiliary routines                                               ====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sort_ja(A)                                                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_matrix), intent(inout) :: A
    ! local variables
    integer :: i,j,k,n
    real(kind(1d0)) :: Aij

    do i=1,A%nrow
        do j=A%ia(i)+1,A%ia(i+1)-1
            n = A%ja(j)
            Aij = A%val(j)
            do k=j-1,A%ia(i),-1
                if (A%ja(k)<=n) exit
                A%ja(k+1) = A%ja(k)
                A%val(k+1) = A%val(k)
            enddo
            A%ja(k+1) = n
            A%val(k+1) = Aij
        enddo
    enddo

end subroutine sort_ja



end module csr_matrix_mod

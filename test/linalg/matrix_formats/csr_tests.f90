program csr_tests

    use csr
    use matrix

    implicit none

    ! sparse matrix
    class(sparse_matrix), allocatable :: A

    ! variables for building the sparse matrix
    integer, allocatable :: rows(:),cols(:)
    real(kind(1d0)), allocatable :: vals(:)

    ! variables for testing accessors/mutators
    integer :: row_indx(3),col_indx(3)
    real(kind(1d0)) :: As(3,3),Bs(3,3)

    ! variables for testing matvec
    real(kind(1d0)), allocatable :: x(:),y(:)

    ! permutation


    ! some other locals
    logical :: match
    integer :: i,j,k
    real(kind(1d0)) :: q



!--------------------------------------------------------------------------!
! Test building the matrix                                                 !
!--------------------------------------------------------------------------!
    allocate(csr_matrix::A)
    call A%init(100,100,298)

    allocate(rows(298),cols(298))
    do i=1,99
        rows(i) = i
        cols(i) = i

        rows(i+100) = i
        cols(i+100) = i+1

        rows(i+199) = i+1
        cols(i+199) = i
    enddo

    rows(100) = 100
    cols(100) = 100

    call A%build(rows,cols)

    do i=1,99
        call A%set_value(i,i,2.d0)
        call A%set_value(i,i+1,-1.d0)
        call A%set_value(i+1,i,-1.d0)
    enddo
    call A%set_value(100,100,2.d0)

    q = A%get_value(100,100)
    if (q/=2.d0) then
        print *, 'A(100,100) should be 2.0'
        print *, 'Value found: ',q
        call exit(1)
    endif

    call A%set_value(1,100,-1.d0)
    call A%set_value(100,1,-1.d0)

    q = A%get_value(100,1)
    if (q/=-1.d0) then
        print *, 'A(100,1) should be -1.0'
        print *, 'Value found: ',q
        call exit(1)
    endif

    row_indx = [ 1,  2, 3 ]
    col_indx = [50, 51, 52]
    As(1,:) = [1.d0, 2.d0, 3.d0]
    As(2,:) = [4.d0, 5.d0, 6.d0]
    As(3,:) = [7.d0, 8.d0, 9.d0]
    call A%set_values(row_indx,col_indx,As)

    Bs = A%get_values(row_indx,col_indx)
    match = .true.
    do j=1,3
        do i=1,3
            if (As(i,j)/=Bs(i,j)) then
                match=.false.
            endif
        enddo
    enddo

    if (.not.match) then
        print *, 'A(1:3,50:52) should be: '
        print *, '  [ 1.0, 2.0, 3.0 ]'
        print *, '  [ 4.0, 5.0, 6.0 ]'
        print *, '  [ 7.0, 8.0, 9.0 ]'
        print *, 'Values found: '
        print *, Bs(1,:)
        print *, Bs(2,:)
        print *, Bs(3,:)
        call exit(1)
    endif



!--------------------------------------------------------------------------!
! Test matrix-vector multiplication                                        !
!--------------------------------------------------------------------------!




!--------------------------------------------------------------------------!
! Test solving triangular systems                                          !
!--------------------------------------------------------------------------!



end program csr_tests

program ellpack_tests

    use ellpack
    use matrix

    implicit none

    ! sparse matrix
    class(sparse_matrix), allocatable :: A

    ! variables for building the sparse matrix
    integer, allocatable :: rows(:),cols(:)
    real(kind(1d0)), allocatable :: vals(:)

    ! variables for testing accessors/mutators
    real(kind(1d0)) :: As(3,3),Bs(3,3)

    ! variables for testing matvec, forward/back solves
    real(kind(1d0)), allocatable :: x(:),y(:)

    ! permutation


    ! some other locals
    integer :: i,j,k,l,row_index(3),col_index(3)
    logical :: match
    real(kind(1d0)) :: q



!--------------------------------------------------------------------------!
! Test building the matrix                                                 !
!--------------------------------------------------------------------------!
    allocate(ellpack_matrix::A)
    call A%init(100,100,300)
    if (A%max_degree/=3) then
        print *, 'Maximum degree of A should be:  ',3
        print *, 'Value found: ',A%max_degree
        call exit(1)
    endif

    allocate(rows(300),cols(300),vals(300))
    vals = 0.d0

    do i=1,100
        rows(i) = i
        cols(i) = i

        rows(i+100) = i
        cols(i+100) = mod(i,100)+1

        rows(i+200) = i
        cols(i+200) = i-1
    enddo
    cols(201) = 100

    call A%build(rows,cols)



!--------------------------------------------------------------------------!
! Test setting and getting entries from the matrix                         !
!--------------------------------------------------------------------------!
    As(1,:) = [  2.d0, -1.d0,  0.d0 ]
    As(2,:) = [ -1.d0,  2.d0, -1.d0 ]
    As(3,:) = [  0.d0, -1.d0,  2.d0 ]
    row_index = [ 1, 2, 3 ]
    col_index = row_index
    call A%set_values(row_index,col_index,As)

    do i=4,100
        call A%set_value(i,i,2.d0)
        call A%set_value(i-1,i,-1.d0)
        call A%set_value(i,i-1,-1.d0)
    enddo

    Bs = A%get_values(row_index,col_index)
    match = .true.
    do j=1,3
        do i=1,3
            if (As(i,j)/=Bs(i,j)) match = .false.
        enddo
    enddo

    if (.not.match) then
        print *, 'A(1:3,1:3) should be: '
        print *, '    [  2.0 -1.0  0.0 ]'
        print *, '    [ -1.0  2.0 -1.0 ]'
        print *, '    [  0.0 -1.0  2.0 ]'
        print *, ' Values found: '
        print *, Bs(1,:)
        print *, Bs(2,:)
        print *, Bs(3,:)

        call exit(1)
    endif

    call A%set_value(1,100,-1.d0)
    call A%set_value(100,1,-1.d0)

    q = A%get_value(1,100)
    if (q/=-1.d0) then
        print *, 'A(1,100) should be -1.0'
        print *, 'Value found: ',q

        call exit(1)
    endif


!--------------------------------------------------------------------------!
! Test matrix multiplication                                               !
!--------------------------------------------------------------------------!
    allocate(x(100),y(100))

    x = 1.d0
    call A%matvec(x,y)

    if (maxval(dabs(y))/=0.d0) then
        print *, 'A*[1.0, ..., 1.0] should be 0.0'
        print *, 'Values found: ', maxval(dabs(y))
    endif

end program ellpack_tests

program bsr_tests

    use bsr
    use matrix

    implicit none

    ! sparse matrix
    class(sparse_matrix), allocatable :: A

    ! variables for building the sparse matrix
    integer, allocatable :: rows(:),cols(:)
    real(kind(1d0)), allocatable :: vals(:)

    ! variables for testing accessors/mutators
    integer :: nbrs(3)
    real(kind(1d0)) :: q,Aij(4,4),Bij(4,4)

    ! vectors for testing matvec, forward/back solves
    real(kind(1d0)), allocatable :: x(:),b(:),z(:)

    ! permutation
    integer, allocatable :: p(:)

    ! some other locals
    integer :: i,j,k,l,row_index(4),col_index(4)
    logical :: match


!--------------------------------------------------------------------------!
! Test building the matrix                                                 !
!--------------------------------------------------------------------------!
    allocate(bsr_matrix::A)
    call A%init(100,100,1200,params=[4,4])

    allocate(rows(75),cols(75))
    vals = 0.d0

    do i=1,25
        rows(i) = i
        cols(i) = i

        rows(i+25) = i
        cols(i+25) = mod(i,25)+1

        rows(i+50) = i
        cols(i+50) = mod(i+23,25)+1
    enddo

    call A%build(rows,cols)

    select type(A)
        type is(bsr_matrix)
            if (A%nbr/=4 .or. A%nbc/=4) then
                print *, 'Rows, cols per block should be:  ',4,4
                print *, 'Values found: ',A%nbr,A%nbc
                call exit(1)
            endif
            if (A%mrow/=25 .or. A%mcol/=25) then
                print *, '# of row, col blocks should be: ',25,25
                print *, 'Values found: ',A%mrow,A%mcol
                call exit(1)
            endif
            if (A%nnz/=1200 .or. A%nnzb/=75) then
                print *, '# of non-zero entries should be: ',1200
                print *, 'Values found: ',A%nnz
                print *, '# of non-zero blocks should be: ',75
                print *, 'Values found: ',A%nnzb
                call exit(1)
            endif
            if (A%max_degree/=3) then
                print *, 'maximum degree should be: ',3
                print *, 'Values found: ',A%max_degree
                call exit(1)
            endif
    end select



!--------------------------------------------------------------------------!
! Test setting and getting entries from the matrix                         !
!--------------------------------------------------------------------------!
    Aij(1,:) = [ 6.d0, -1.d0, -1.d0,  0.d0 ]
    Aij(2,:) = [-1.d0,  6.d0, -1.d0, -1.d0 ]
    Aij(3,:) = [-1.d0, -1.d0,  6.d0, -1.d0 ]
    Aij(4,:) = [ 0.d0, -1.d0, -1.d0,  6.d0 ]
    Bij = Aij
    do i=1,25
        j = 4*(i-1)+1
        row_index = [ j, j+1, j+2, j+3 ]
        col_index = row_index
        call A%set_values(row_index,col_index,Aij)
    enddo

    Aij(1,:) = [-1.d0,  0.d0,  0.d0,  0.d0 ]
    Aij(2,:) = [ 0.d0, -1.d0,  0.d0,  0.d0 ]
    Aij(3,:) = [-1.d0,  0.d0, -1.d0,  0.d0 ]
    Aij(4,:) = [-1.d0, -1.d0,  0.d0, -1.d0 ]
    do i=1,25
        j = 4*(i-1)+1
        k = mod( 4*i, 100 )+1
        row_index = [ j, j+1, j+2, j+3 ]
        col_index = [ k, k+1, k+2, k+3 ]
        call A%set_values(row_index,col_index,Aij)
        call A%set_values(col_index,row_index,transpose(Aij))
    enddo

    q = A%get_value(1,1)
    if (q/=6.d0) then
        print *, 'A(1,1) should be: ',6.d0
        print *, 'Values found: ',q
        call exit(1)
    endif

    q = A%get_value(1,2)
    if (q/=-1.d0) then
        print *, 'A(1,2) should be: ',-1.d0
        print *, 'Values found: ',q
        call exit(1)
    endif

    Aij = A%get_values([1,2,3,4],[9,10,11,12])
    if (maxval(dabs(Aij))/=0.d0) then
        print *, 'A(1:4,9:12) should be 0'
        print *, 'Values found: '
        print *, Aij(1,:)
        print *, Aij(2,:)
        print *, Aij(3,:)
        print *, Aij(4,:)
        call exit(1)
    endif


    Aij = A%get_values([1,2,3,4],[1,2,3,4])

    Bij(1,:) = [ 6.d0, -1.d0, -1.d0,  0.d0 ]
    Bij(2,:) = [-1.d0,  6.d0, -1.d0, -1.d0 ]
    Bij(3,:) = [-1.d0, -1.d0,  6.d0, -1.d0 ]
    Bij(4,:) = [ 0.d0, -1.d0, -1.d0,  6.d0 ]

    match = .true.
    do j=1,4
        do i=1,4
            if (Aij(i,j)/=Bij(i,j)) match = .false.
        enddo
    enddo

    if (.not.match) then
        print *, 'A(1:4,1:4) should be: '
        print *, '    [  6.0 -1.0 -1.0  0.0  ]'
        print *, '    [ -1.0  6.0 -1.0 -1.0  ]'
        print *, '    [ -1.0 -1.0  6.0 -1.0  ]'
        print *, '    [  0.0 -1.0 -1.0  6.0  ]'
        print *, 'Values found: '
        print *, Aij(1,:)
        print *, Aij(2,:)
        print *, Aij(3,:)
        print *, Aij(4,:)

        call exit(1)
    endif


    Aij = A%get_values([5,6,7,8],[9,10,11,12])

    Bij(1,:) = [-1.d0,  0.d0,  0.d0,  0.d0 ]
    Bij(2,:) = [ 0.d0, -1.d0,  0.d0,  0.d0 ]
    Bij(3,:) = [-1.d0,  0.d0, -1.d0,  0.d0 ]
    Bij(4,:) = [-1.d0, -1.d0,  0.d0, -1.d0 ]

    match = .true.
    do j=1,4
        do i=1,4
            if (Aij(i,j)/=Bij(i,j)) match = .false.
        enddo
    enddo

    if (.not.match) then
        print *, 'A(5:8,9:12) should be: '
        print *, '    [ -1.0  0.0  0.0  0.0  ]'
        print *, '    [  0.0 -1.0  0.0  0.0  ]'
        print *, '    [ -1.0  0.0 -1.0  0.0  ]'
        print *, '    [ -1.0 -1.0  0.0 -1.0  ]'
        print *, 'Values found: '
        print *, Aij(1,:)
        print *, Aij(2,:)
        print *, Aij(3,:)
        print *, Aij(4,:)

        call exit(1)
    endif



!--------------------------------------------------------------------------!
! Test matrix multiplication                                               !
!--------------------------------------------------------------------------!
    allocate(x(100),z(100))
    x = 1.d0
    z = 1.d0
    call A%matvec(x,z)

    if (maxval(dabs(z))/=0) then
        print *, 'A*[1, 1, ... , 1] should be 0'
        print *, 'Values found: ', maxval(dabs(z))
    endif



!--------------------------------------------------------------------------!
! Test solution of triangular systems                                      !
!--------------------------------------------------------------------------!



!--------------------------------------------------------------------------!
! Test converting to coordinate format                                     !
!--------------------------------------------------------------------------!
    deallocate(rows,cols)
    allocate(rows(1200),cols(1200),vals(1200))
    call A%convert_to_coo(rows,cols,vals)

end program bsr_tests

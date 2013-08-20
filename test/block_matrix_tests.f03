program block_matrix_tests

use fempack

implicit none

    class(graph), allocatable :: g
    class(block_sparse_matrix), allocatable :: A
    integer :: i,j,k
    integer, allocatable :: edges(:,:), nbrs(:)
    real(dp), allocatable :: block(:,:), x(:,:), y(:,:), z(:)

    allocate(edges(2,31))

    edges(:,1)  = [1, 2]
    edges(:,2)  = [1, 3]
    edges(:,3)  = [1, 4]
    edges(:,4)  = [1, 5]
    edges(:,5)  = [1, 6]
    edges(:,6)  = [1, 7]
    edges(:,7)  = [2, 3]
    edges(:,8)  = [3, 4]
    edges(:,9)  = [4, 5]
    edges(:,10) = [5, 6]
    edges(:,11) = [6, 7]
    edges(:,12) = [7, 2]

    do i=1,12
        edges(1,i+12) = edges(2,i)
        edges(2,i+12) = edges(1,i)
    enddo

    do i=1,7
        edges(:,i+24) = [i, i]
    enddo

    allocate(cs_graph::g)
    allocate(bcsr_matrix::A)

    call g%init(7,7,edges)
    call A%init(14,14)
    call A%assemble(g)

    allocate(block(2,2), nbrs(8))

    block(:,1) = [ 2.0_dp, -1.0_dp]
    block(:,2) = [-1.0_dp,  2.0_dp]

    do i=1,7
        call A%neighbors(i,nbrs)
        do k=1,A%max_degree
            j = nbrs(k)
            if (j/=0 .and. j/=i) then
                call A%set_block(i,j,-block)
                call A%add_block(i,i,block)
            endif
        enddo
    enddo

    call A%get_block(1,1,block)
    if ( .not.(block(1,1)==12.0_dp .and. block(2,1)==-6.0_dp &
        & .and. block(1,2)==-6.0_dp .and. block(2,2)==12.0_dp) ) then
        print *, 'A(1,1) should be = [12.0, -6.0; -6.0, 12.0]'
        print *, 'Values found:  ',block(1,:)
        print *, '               ',block(2,:)
    endif

    call A%get_block(1,2,block)
    if ( .not.(block(1,1)==-2.0_dp .and. block(2,1)==1.0_dp &
        & .and. block(1,2)==1.0_dp .and. block(2,2)==-2.0_dp) ) then
        print *, 'A(1,2) should be = [-2.0, 1.0; 1.0, -2.0]'
        print *, 'Values found:  ',block(1,:)
        print *, '               ',block(2,:)
    endif

    allocate(x(2,7),y(2,7),z(14))
    x = 1.0_dp
    y = 1.0_dp
    call A%matmul(x,y)
    call A%matmul(x,z)

    if ( maxval(dabs(y))>1.0e-14 ) then
        print *, 'A%[1,...,1] should be = 0;'
        print *, 'max(abs(y)) = ',maxval(dabs(y))
    endif

    if ( maxval(dabs(z))>1.0e-14 ) then
        print *, 'A%[1,...,1] should be = 0;'
        print *, 'max(abs(z)) = ',maxval(dabs(z))
    endif

    deallocate(g,A,x,y,z,edges)


end program block_matrix_tests

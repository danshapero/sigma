program matrix_tests

use fempack

implicit none

    class(graph), allocatable :: g
    class(sparse_matrix), allocatable :: A
    integer :: i,j,k,degree
    real(dp) :: z
    integer, allocatable :: edges(:,:), nbrs(:)
    real(dp), allocatable :: x(:), y(:)

    allocate(edges(2,24))

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

    allocate(cs_graph::g)
    call g%init(7,7,edges)

    allocate(csr_matrix::A)
    call A%assemble(g)

    allocate(nbrs(7))

    ! Fill A to be the graph Laplacian
    do i=1,7
        call A%neighbors(i,nbrs)
        degree = 0
        do k=1,A%max_degree
            if (nbrs(k)/=0 .and. nbrs(k)/=i) degree = degree+1
        enddo
        do k=1,A%max_degree
            j = nbrs(k)
            if (j/=0) then
                call A%set_value(i,j,-1.0_dp)
                call A%add_value(i,i,1.0_dp)
            endif
        enddo
    enddo

    z = A%get_value(1,1)
    if (z/=6.0_dp) then
        print *, 'A(1,1) should be = 6.0 for A the graph Laplacian; '
        print *, 'value found: ',z
    endif

    allocate(x(7),y(7))

    x = 1.0_dp
    y = 1.0_dp
    call A%matvec(x,y)

    if ( maxval(dabs(y))>1.0e-14 ) then
        print *, 'A*[1,...,1] should be = 0;'
        print *, 'max(abs(y)) = ',maxval(dabs(y))
    endif

end program matrix_tests

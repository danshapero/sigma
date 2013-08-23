program matrix_tests

use fempack

implicit none

    class(graph), allocatable :: g, h
    class(sparse_matrix), allocatable :: A, B
    integer :: i,j,k,degree,test
    real(dp) :: z
    logical :: correct
    integer, allocatable :: edges(:,:), nbrs(:), p(:)
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

    allocate(nbrs(8),p(7))
    allocate(x(7),y(7))

    do i=1,7
        p(i) = i-1
    enddo
    p(1) = 7

    do test=1,3
        ! Allocate the graph & matrix
        select case(test)
            case(1)
                call new_graph(g,'cs',7,7,edges)
                call new_sparse_matrix(A,'csr',7,7,g)

                call new_graph(h,'cs',7,7)
                call new_sparse_matrix(B,'csr',7,7,h)
            case(2)
                call new_graph(g,'cs',7,7,edges)
                call new_sparse_matrix(A,'csc',7,7,g)

                call new_graph(h,'cs',7,7)
                call new_sparse_matrix(B,'csc',7,7,h)
            case(3)
                call new_graph(g,'coo',7,7,edges)
                call new_sparse_matrix(A,'coo',7,7,g)

                call new_graph(h,'coo',7,7)
                call new_sparse_matrix(B,'coo',7,7,h)
        end select

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
            call B%set_value(i,i,1.0_dp)
        enddo

        ! Test matrix multiplications
        z = A%get_value(1,1)
        if (z/=6.0_dp) then
            print *, 'A(1,1) should be = 6.0 for A the graph Laplacian; '
            print *, 'value found: ',z
        endif

        x = 1.0_dp
        y = 1.0_dp
        call A%matvec(x,y)

        if ( maxval(dabs(y))>1.0e-14 ) then
            print *, 'A*[1,...,1] should be = 0;'
            print *, 'range(y) = ',minval(y),maxval(y)
        endif

        ! Test adding two matrices
        call A%sub_matrix_add(B)
        y = 0.0_dp
        call A%matvec(x,y)
        if ( maxval(dabs(y)-1.0)>1.0e-14 ) then
            print *, '(A+I)*[1,...,1] should be = 1;'
            print *, 'range(y) = ',minval(y),maxval(y)
        endif

        ! Test permuting the matrices
        call A%right_permute(p)
        call A%left_permute(p)
        call A%neighbors(7,nbrs)
        nbrs(1:7) = nbrs(order(nbrs(1:7)))
        correct = .true.
        do i=1,7
            if (nbrs(i)/=i) correct = .false.
        enddo
        if (.not.correct) then
            print *, 'Permutation failed'
        endif

        ! Free up the graphs and matrices for the next test
        deallocate(g,A,h,B)
    enddo

    deallocate(x,y,edges,nbrs)

end program matrix_tests

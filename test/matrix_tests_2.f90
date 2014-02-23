program matrix_tests_2

use sigma

implicit none

    type(sparse_matrix) :: A, B
    class(graph), pointer :: g, h, gr, hr

    integer :: i,j,k,nn,test1,test2
    integer, allocatable :: neighbors(:)

    real(dp) :: p, error
    real(dp), allocatable :: u(:), x(:), y(:), z(:)

    call init_seed()
    nn = 64
    p = 8.0/nn

    !----------------------------------------------------------------------!
    ! Construct reference graphs from which all test graphs are copied     !
    !----------------------------------------------------------------------!
    allocate(ll_graph::gr)
    allocate(ll_graph::hr)
    call gr%init(nn)
    call hr%init(nn)

    allocate(u(nn),x(nn),y(nn),z(nn))
    ! Randomly generate graphs
    do i=1,nn
        call gr%add_edge(i,i)

        call random_number(z)
        call random_number(y)
        ! Connect each pair (i,j) in gr with probability p
        do j=i+1,nn
            if (z(j)>p) then
                call gr%add_edge(i,j)
                call gr%add_edge(j,i)

                ! Connect them in hr with probability p/2
                if (y(j)>0.5_dp) then
                    call hr%add_edge(i,j)
                    call hr%add_edge(j,i)
                endif
            endif
        enddo
    enddo

    allocate(neighbors(gr%max_degree))
    

    !----------------------------------------------------------------------!
    ! Test adding two matrices with each possible connectivity graph       !
    !----------------------------------------------------------------------!
    do test1=1,4
        ! Allocate g to each possible graph type
        select case(test1)
            case(1)
                allocate(ll_graph::g)
            case(2)
                allocate(coo_graph::g)
            case(3)
                allocate(cs_graph::g)
            case(4)
                allocate(ellpack_graph::g)
        end select
        call g%init(gr)

        ! Make a matrix from g
        call A%init(nn,nn,'row',g)

        do test2=1,4
            ! Allocate h to each possible graph type
            select case(test1)
                case(1)
                    allocate(ll_graph::h)
                case(2)
                    allocate(coo_graph::h)
                case(3)
                    allocate(cs_graph::h)
                case(4)
                    allocate(ellpack_graph::h)
            end select
            call h%init(hr)

            ! Make a matrix from h
            call B%init(nn,nn,'row',h)

            ! Make A the graph Laplacian
            call A%zero()
            do i=1,nn
                call A%g%neighbors(i,neighbors)

                do k=1,A%max_degree
                    j = neighbors(k)
                    if (j/=0 .and. j/=i) then
                        call A%set_value(i,j,-1.0_dp)
                        call A%add_value(i,i,1.0_dp)
                    endif
                enddo
            enddo

            ! Make B anti-symmetric
            call B%zero()
            do i=1,nn
                call B%g%neighbors(i,neighbors)

                do k=1,B%max_degree
                    j = neighbors(k)
                    if (j/=0) then
                        if (j>i) call B%set_value(i,j,1.0_dp)
                        if (i<j) call B%set_value(i,j,-1.0_dp)
                    endif
                enddo
            enddo

            ! Make x a bunch of random numbers
            call random_number(x)

            ! Compute y = A*x
            call A%matvec(x,y)

            ! Compute z = B*x
            call B%matvec(x,z)
    
            ! Add B to A
            call A%add(B)

            ! Compute u = A*x
            call A%matvec(x,u)

            error = maxval(dabs(y+z-u))
            if (error>1.0e-12) then
                print *, test1,test2,error
                call exit(1)
            endif

            deallocate(h)
            call B%destroy()
        enddo

        deallocate(g)
        call A%destroy()
    enddo


end program matrix_tests_2

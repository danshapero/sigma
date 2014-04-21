program matrix_tests_4

use sigma

implicit none

    type(sparse_matrix) :: A, B, C
    class(graph), pointer :: g, h, gr, hr

    integer :: i, j, k, nn, test1, test2
    integer, allocatable :: neighbors(:)

    real(dp) :: p, q, error
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
    do i=1,nn
        call random_number(z)
        call random_number(y)

        do j=1,nn
            if (z(j)<p) call gr%add_edge(i,j)
            if (y(j)<p) call hr%add_edge(i,j)
        enddo
    enddo

    allocate(neighbors(max(gr%max_degree,hr%max_degree)))


    do test1=1,4
        print *, test1

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

        ! Make B a random matrix
        call B%init(nn,nn,'row',g)
        do i=1,nn
            call B%g%neighbors(i,neighbors)
            do k=1,B%g%max_degree
                j = neighbors(k)
                if (j/=0) then
                    call random_number(q)
                    call B%set_value(i,j,2*q-1)
                endif
            enddo
        enddo

        do test2=1,4
            print *, test2
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

            ! Make C a random matrix
            call C%init(nn,nn,'row',h)
            do i=1,nn
                call C%g%neighbors(i,neighbors)
                do k=1,C%g%max_degree
                    j = neighbors(k)
                    if (j/=0) then
                        call random_number(q)
                        call C%set_value(i,j,2*q-1)
                    endif
                enddo
            enddo

            ! Add B and C into the matrix A
            call A%add(B,C)

            x = 0.0_dp
            u = 0.0_dp
            y = 0.0_dp
            z = 0.0_dp

            ! Make x a bunch of random numbers
            call random_number(x)

            ! Compute u = A*x
            call A%matvec(x,u)

            ! Compute y = B*x, z = C*x
            call B%matvec(x,y)
            call C%matvec(x,z)

            error = maxval(dabs(y+z-u))
            if (error>1.0e-12) then
                print *, test1,test2,error
            endif

            call C%destroy()
            call A%destroy()
            deallocate(h)
        enddo

        call B%destroy()
        deallocate(g)
    enddo

end program matrix_tests_4

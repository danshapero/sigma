!--------------------------------------------------------------------------!
program matrix_tests_4                                                     !
!--------------------------------------------------------------------------!
!     This program tests creating a new matrix as the sum of two other     !
! random matrices with each possible graph type.                           !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! Matrices and graphs
    type(sparse_matrix) :: A, B, C
    class(graph), pointer :: g, h, gr, hr, gph
    character(len=3) :: orientation1, orientation2, orientation3
    ! Graph edge iterators
    type(graph_edge_cursor) :: cursor
    integer :: num_batches, num_returned, edges(2,batch_size)
    ! Integer indices
    integer :: i, j, k, n, nn, test1, test2, test3, frmt1, frmt2, frmt3
    ! Random numbers and vectors
    real(dp) :: p, q, error
    real(dp), allocatable :: u(:), x(:), y(:), z(:)
    ! command-line arguments
    character(len=16) :: arg
    logical verbose
    ! strings for printing verbose stuff
    character(len=24) :: message


    ! Get command line arguments to see if we're running in verbose mode
    verbose = .false.
    call getarg(1,arg)
    select case(trim(arg))
        case("-v")
            verbose = .true.
        case("-V")
            verbose = .true.
        case("--verbose")
            verbose = .true.
    end select


    ! Initialize a random seed
    call init_seed()
    nn = 72
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


    !----------------------------------------------------------------------!
    ! Outer loop: construct a random matrix of each type & order           !
    !----------------------------------------------------------------------!
    do test1=1,4
    do frmt1=1,2
        ! Select the orientation of the matrix
        select case(frmt1)
            case(1)
                orientation1 = "row"
            case(2)
                orientation1 = "col"
        end select

        ! Allocate g to each possible graph type
        select case(test1)
            case(1)
                allocate(ll_graph::g)
                message = "linked-list"
            case(2)
                allocate(coo_graph::g)
                message = "coordinate"
            case(3)
                allocate(cs_graph::g)
                message = "compressed sparse"
            case(4)
                allocate(ellpack_graph::g)
                message = "ellpack"
        end select

        if (verbose) then
            print *, "Testing ",trim(message)," graph with ", &
                & orientation1," orientation"
        endif

        call g%init(gr)

        ! Make B a random matrix
        call B%init(nn,nn,orientation1,g)
        cursor = B%g%make_cursor(0)
        num_batches = (cursor%start-cursor%final)/batch_size+1
        do n=1,num_batches
            call B%g%get_edges(edges,cursor,batch_size,num_returned)

            do k=1,num_returned
                i = edges(B%order(1),k)
                j = edges(B%order(2),k)

                if (i/=0 .and. i/=0) then
                    call random_number(q)
                    call B%set_value(i,j,2*q-1)
                endif
            enddo
        enddo


        !------------------------------------------------------------------!
        ! Inner loop: construct another random matrix of each type & order !
        !------------------------------------------------------------------!
        do test2=1,4
        do frmt2=1,2
            ! Select the orientation of the matrix
            select case(frmt2)
                case(1)
                    orientation2 = "row"
                case(2)
                    orientation2 = "col"
            end select

            ! Allocate h to each possible graph type
            select case(test2)
                case(1)
                    allocate(ll_graph::h)
                    message = "linked-list"
                case(2)
                    allocate(coo_graph::h)
                    message = "coordinate"
                case(3)
                    allocate(cs_graph::h)
                    message = "compressed sparse"
                case(4)
                    allocate(ellpack_graph::h)
                    message = "ellpack"
            end select

            if (verbose) then
                print *, "  Testing ",trim(message), " graph with ", &
                    & orientation2," orientation"
            endif

            call h%init(hr)

            ! Make C a random matrix
            call C%init(nn,nn,orientation2,h)
            cursor = C%g%make_cursor(0)
            num_batches = (cursor%final-cursor%start)/batch_size+1
            do n=1,num_batches
                call C%g%get_edges(edges,cursor,batch_size,num_returned)

                do k=1,num_returned
                    i = edges(C%order(1),k)
                    j = edges(C%order(2),k)

                    call random_number(q)
                    call C%set_value(i,j,2*q-1)
                enddo
            enddo


            !--------------------------------------------------------------!
            ! Check that the sum of these matrices is computed correctly   !
            !--------------------------------------------------------------!
            do test3=1,4
            do frmt3=1,2
                ! Select the orientation of the matrix
                select case(frmt3)
                    case(1)
                        orientation3 = "row"
                    case(2)
                        orientation3 = "col"
                end select

                ! Allocate gph to each possible graph type
                select case(test3)
                    case(1)
                        allocate(ll_graph::gph)
                        message = "linked-list"
                    case(2)
                        allocate(coo_graph::gph)
                        message = "coordinate"
                    case(3)
                        allocate(cs_graph::gph)
                        message = "compressed sparse"
                    case(4)
                        allocate(ellpack_graph::gph)
                        message = "ellpack"
                end select

                if (verbose) then
                    print *, "    Testing ",trim(message), " graph with ", &
                        & orientation3," orientation"
                endif


                ! Add B and C into the matrix A
                call add_sparse_matrices(A, B, C, &
                    & g = gph, orientation = orientation3)

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
                    print *, 'On test',test1,test2,test3
                    print *, 'Matrix sum A = B+C failed; should have'
                    print *, 'A*x = B*x+C*x, but error is',error
                    print *, 'Terminating.'
                    call exit(1)
                endif

                call A%destroy()
                call gph%destroy()
                deallocate(gph)
            enddo
            enddo

            call C%destroy()
            call h%destroy()
            deallocate(h)
        enddo
        enddo

        call B%destroy()
        call g%destroy()
        deallocate(g)
    enddo
    enddo

end program matrix_tests_4

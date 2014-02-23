program graph_tests_2

use sigma

implicit none

    class(graph), allocatable :: g, h
    integer :: i, j, k, l, n, degree, test, min_degree
    real(dp) :: z(16)
    ! variables for graph edge iteration
    integer :: num_blocks, num_returned, edges(2,64)
    type(graph_edge_cursor) :: cursor
    ! command-line arguments
    character(len=16) :: arg
    logical :: verbose



    ! Get command line arguments to see if we're running in verbose mode
    verbose = .false.
    call getarg(1,arg)
    if (trim(arg)=="-v" .or. trim(arg)=="--verbose") then
        verbose = .true.
    endif

    ! Initialize a random seed
    call init_seed()


    !----------------------------------------------------------------------!
    ! Create a random sparse graph                                         !
    !----------------------------------------------------------------------!
    allocate(ll_graph::h)
    call h%init(512,512)

    do i=1,512
        ! Initialize the degree of node i to be 0
        degree = 0

        ! Make some random numbers
        call random_number(z)

        ! degree = number of samples greater than 1/2
        do k=1,16
            if (z(k)>0.5_dp) degree = degree+1
        enddo

        ! Make some more random numbers
        call random_number(z)

        do k=1,degree
            ! Get a random number j between 1 and 512
            j = min(int(z(k)*512)+1,512)

            ! Make (i,j) connected in h
            call h%add_edge(i,j)
            call h%add_edge(j,i)
        enddo
    enddo

    min_degree = h%max_degree

    select type(h)
        type is(ll_graph)
            do i=1,512
                min_degree = min( h%lists(i)%length, min_degree )
            enddo
    end select

    if (verbose) then
        print *, 'Random graph generated'
        print *, 'Min/max degree: ',min_degree,h%max_degree
        print *, ' '
    endif


    !----------------------------------------------------------------------!
    ! Copy the reference graph h to graphs in other formats                !
    !----------------------------------------------------------------------!
    if (verbose) print *, 'Testing copy constructor'

    do test=1,4
        select case(test)
            case(1)
                allocate(ll_graph::g)
                if (verbose) print *, '  Test 1: linked-list graph'
            case(2)
                allocate(coo_graph::g)
                if (verbose) print *, '  Test 2: coordinate graph'
            case(3)
                allocate(cs_graph::g)
                if (verbose) print *, '  Test 3: compressed sparse graph'
            case(4)
                allocate(ellpack_graph::g)
                if (verbose) print *, '  Test 4: ellpack graph'
        end select

        call g%copy(h)

        ! Iterate through all the edges of g and make sure they all are
        ! connected in h
        cursor = g%make_cursor(0)
        num_blocks = (cursor%final-cursor%start)/64+1

        do n=1,num_blocks
            edges = g%get_edges(cursor,64,num_returned)

            do k=1,num_returned
                i = edges(1,k)
                j = edges(2,k)

                if (i/=0 .and. j/=0) then
                    l = h%find_edge(i,j)
                    if (l==-1) then
                        print *, '    On test',test,'vertices',i,j
                        print *, '    not connected in reference graph,'
                        print *, '    but are connected in copied graph.'
                        call exit(1)
                    endif
                endif
            enddo
        enddo

        ! Iterate through all the edges of h and make sure they all are
        ! connected in g
        cursor = h%make_cursor(0)
        num_blocks = (cursor%final-cursor%start)/64+1

        do n=1,num_blocks
            edges = h%get_edges(cursor,64,num_returned)

            do k=1,num_returned
                i = edges(1,k)
                j = edges(2,k)

                if (i/=0 .and. j/=0) then
                    l = g%find_edge(i,j)
                    if (l==-1) then
                        print *, '    On test',test,'vertices',i,j
                        print *, '    are connected in reference graph,'
                        print *, '    but not connected in copied graph.'
                        call exit(1)
                    endif
                endif
            enddo
        enddo

        deallocate(g)
    enddo



    !----------------------------------------------------------------------!
    ! Make another random sparse graph and add the two graphs              !
    !----------------------------------------------------------------------!


end program graph_tests_2

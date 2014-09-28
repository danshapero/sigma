!--------------------------------------------------------------------------!
program graph_test_basics                                                  !
!--------------------------------------------------------------------------!
! This program performs tests of basic graph operations:                   !
!     o initialization                                                     !
!     o adding and removing edges                                          !
!     o checking the connectedness of two vertices                         !
!     o returning all neighbors of a vertex                                !
!     o edge iterators                                                     !
!     o permutation                                                        !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph
    class(graph_interface), pointer :: g

    ! dense graph
    integer, allocatable :: A(:,:), B(:,:), BP(:,:)

    ! integer indices
    integer :: i, j, k, max_degree, nn, ne, frmt

    ! random numbers
    real(dp) :: z, c

    ! variables for getting neighbors
    integer :: degree
    integer, allocatable :: neighbors(:)

    ! variables for iterating over graph edges
    integer :: n, num_batches, num_returned, edges(2, batch_size)
    type(graph_edge_cursor) :: cursor

    ! permutation
    integer, allocatable :: p(:)

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose

    ! miscellaneous
    logical :: con



    !----------------------------------------------------------------------!
    ! Get command line arguments to see if we're running in verbose mode   !
    !----------------------------------------------------------------------!
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



    !----------------------------------------------------------------------!
    ! Set the graph size and initialize a random seed                      !
    !----------------------------------------------------------------------!

    nn = 64
    c = log(1.0_dp * nn) / log(2.0_dp) / nn

    call init_seed()



    !----------------------------------------------------------------------!
    ! Make a random reference graph, stored in a dense array               !
    !----------------------------------------------------------------------!

    allocate(A(nn, nn), B(nn, nn))
    A = 0
    B = 0

    do j = 1, nn
        do i = 1, nn
            call random_number(z)
            if (z < c) A(i, j) = 1
        enddo
    enddo

    ! Compute the max degree of the graph
    max_degree = 0
    ne = 0
    do i = 1, nn
        k = count(A(i, :) /= 0)

        if (k == 0) then
            call random_number(z)
            j = int(z * nn) + 1
            A(i, j) = 1
            k = 1
        endif

        ne = ne + k
        max_degree = max(max_degree, k)
    enddo

    allocate(neighbors(max_degree))



    !----------------------------------------------------------------------!
    ! Make a random permutation and permute the reference graph            !
    !----------------------------------------------------------------------!

    allocate(p(nn))
    do i = 1, nn
        p(i) = i
    enddo

    do i = nn, 2, -1
        ! Pick a random number `j` between 1 and `i`
        call random_number(z)
        j = int(z * i) + 1

        ! Swap `p(i)` and `p(j)`
        k = p(i)
        p(i) = p(j)
        p(j) = k
    enddo

    allocate(BP(nn, nn))



    !----------------------------------------------------------------------!
    ! Test each graph type                                                 !
    !----------------------------------------------------------------------!
    do frmt = 1, num_graph_types
        B = 0

        if (verbose) print *, 'Format #', frmt

        ! Allocate the graph
        call choose_graph_type(g, frmt)

        ! Initialize the graph, knowing that the maximum degree is `d`
        call g%init(nn, nn)


        ! Check that the graph has the right number of nodes and edges
        if (g%n /= nn) then
            print *, '    Initializing graph failed, number of nodes'
            print *, '    should be:   ', nn
            print *, '    Number found:', g%n
        endif

        if (g%ne /= 0) then
            print *, '    Just-initialized graph should have 0 edges'
            print *, '    Number of edges found:', g%ne
            call exit(1)
        endif


        ! Add in the edges from the dense array
        do j = 1, nn
            do i = 1, nn
                if (A(i, j) /= 0) then
                    call g%add_edge(i, j)
                endif
            enddo
        enddo


        ! Check that `g` has the right number of edges now
        if (g%ne /= ne) then
            print *, '    After inserting edges, total number should be', ne
            print *, '    Number of edges found:', g%ne
            call exit(1)
        endif


        ! Check that the `connected` method works
        do j = 1, nn
            do i = 1, nn
                con = (A(i,j) == 1)
                if (g%connected(i, j) .neqv. con) then
                    print *, '    Graph returned incorrect connectivity for'
                    print *, '    edges:', i, j
                    print *, '    Connected in g:', g%connected(i, j)
                    print *, '    Connected in A:', con
                    call exit(1)
                endif
            enddo
        enddo


        ! Check that the graph's maximum degree is correct
        if (g%max_degree() /= max_degree) then
            print *, '    Max degree of g should be:', max_degree
            print *, '    Max degree found:', g%max_degree()
            call exit(1)
        endif


        ! Check that converting the graph to a dense array works properly
        call g%to_dense_graph(B)
        if ( maxval(abs(A - B)) /= 0 ) then
            print *, '    Converting graph to dense array failed.'
            call exit(1)
        endif

        call g%to_dense_graph(B, trans = .true.)
        if ( maxval(abs(transpose(A) - B)) /= 0 ) then
            print *, '    Converting graph transpose to dense array failed.'
            call exit(1)
        endif

        call g%to_dense_graph(B)


        ! Check that finding the degree of a vertex works
        do i = 1, nn
            degree = count( A(i, :) /= 0 )
            if (g%degree(i) /= degree) then
                print *, '    Getting degree failed.'
                print *, '    Degree of vertex', i ,'should be:', degree
                print *, '    Degree found:', g%degree(i)
                call exit(1)
            endif
        enddo


        ! Check that finding all neighbors of a vertex works
        do i = 1, nn
            ! Get all the neighbors of vertex `i`
            call g%get_neighbors(neighbors, i)

            do k = 1, max_degree
                ! For each neighbor `j` of `i`,
                j = neighbors(k)

                ! Toggle the value of B(i, j).
                if (j /= 0) B(i, j) = A(i, j) - 1
            enddo
        enddo

        ! If the neighbors were returned right, `B` should be zero.
        ! If the graph returned a vertex `j` of some vertex `i` which is
        ! not actually a neighbor of `i`, then B(i, j) = -1.
        ! Conversely, if a neighbor `j` of `i` was not returned, then
        ! B(i, j) = 1.
        if ( maxval(abs(B)) /= 0 ) then
            print *, '    Getting neighbors failed.'

            if (maxval(B) > 0) then
                print *, '    Not all neighbors were returned.'
            endif

            if (minval(B) < 0) then
                print *, '    False positive neighbors were returned.'
            endif

            call exit(1)
        endif

        ! Convert g to a dense graph again
        call g%to_dense_graph(B)


        ! Check that iterating through the graph's edges works
        cursor = g%make_cursor()
        num_batches = (cursor%last - cursor%start) / batch_size + 1

        do n = 1, num_batches
            call g%get_edges(edges, cursor, batch_size, num_returned)

            do k = 1, num_returned
                i = edges(1, k)
                j = edges(2, k)

                ! Check that the edge returned is not null
                if (i == 0 .or. j == 0) then
                    print *, '    Graph edge iterator failed; edge #', k
                    print *, '    returned a 0 vertex.'
                    call exit(1)
                endif

                B(i, j) = A(i, j) - 1
            enddo
        enddo


        ! If the edge iterator worked properly, `B` should be zero;
        ! see comment above for neighbors check for explanation.
        if ( maxval(abs(B)) /= 0 ) then
            print *, '    Iterating over all graph edges failed.'

            if (maxval(B) > 0) then
                print *, '    Not all edges were returned.'
            endif

            if (minval(B) < 0) then
                print *, '    False positive edges were returned.'
            endif

            call exit(1)
        endif


        ! Check that permutations work
        call g%left_permute(p)
        call g%right_permute(p)

        do j = 1, nn
            do i = 1, nn
                BP(p(i), p(j)) = A(i, j)
            enddo
        enddo

        call g%to_dense_graph(B)

        if ( maxval(abs(B - BP)) > 0 ) then
            print *, '    Graph permutation failed.'
            call exit(1)
        endif


        ! Check that deleting edges works.
        do i = 1, nn
            ! First, find a vertex with degree greater than 2,
            degree = g%degree(i)

            if (degree > 2) then
                ! then delete one of its neighbors
                call g%get_neighbors(neighbors, i)

                j = neighbors(1)
                call g%delete_edge(i, j)
                BP(i, j) = 0
                exit
            endif
        enddo

        call g%to_dense_graph(B)

        if ( maxval(abs(B - BP)) > 0 ) then
            print *, '    Deleting edges failed.'
            call exit(1)
        endif


        ! Clear the graph data structure for the next test.
        call g%destroy()
        deallocate(g)

    enddo   ! End of loop over `frmt`


    deallocate(A, B, BP, p, neighbors)

end program graph_test_basics

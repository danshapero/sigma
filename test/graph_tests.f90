program graph_tests

use fempack

implicit none

    ! variables for testing correctness of graph operations
    class(graph), allocatable :: g
    integer :: i,j,k,l,test
    integer, allocatable :: edges(:,:), nbrs(:), p(:)
    logical :: correct
    ! variables for testing graph edge iterator
    type(graph_edge_cursor) :: cursor
    integer, allocatable :: edg(:,:)
    integer :: num_returned
    logical, allocatable :: found_by_iterator(:)

    allocate(edges(2,12))

    edges(:,1) = [1, 2]
    edges(:,2) = [1, 3]
    edges(:,3) = [1, 4]
    edges(:,4) = [1, 5]
    edges(:,5) = [1, 6]
    edges(:,6) = [1, 7]
    edges(:,7) = [2, 3]
    edges(:,8) = [3, 4]
    edges(:,9) = [4, 5]
    edges(:,10) = [5, 6]
    edges(:,11) = [6, 7]
    edges(:,12) = [7, 2]

    allocate(edg(2,12), found_by_iterator(24))


    !----------------------------------------------------------------------!
    ! Loop through and test each graph type                                !
    !----------------------------------------------------------------------!
    do test=1,4
        ! Allocate the graph
        select case(test)
            case(1)
                allocate(ll_graph::g)
            case(2)
                allocate(coo_graph::g)
            case(3)
                allocate(cs_graph::g)
            case(4)
                allocate(ellpack_graph::g)
        end select

        ! Initialize the graph with a set of pre-defined edges
        call g%init(7,7,edges)

        ! Check that the initialization routine has correctly entered the
        ! number of edges in the graph
        if (g%ne/=12) then
            print *, 'Graph should have 12 edges'
            call exit(1)
        endif

        ! Check that edge 1 has been connected to edge 2 but that edge 2
        ! is not connected to edge 1
        if (.not.g%connected(1,2) .or. g%connected(2,1)) then
            print *, 'Graph should have 1 -> 2, 2 -/> 1; we have:'
            if (.not.g%connected(1,2)) then
                print *, ' 1 -/> 2'
            else
                print *, ' 2 -> 1'
            endif
            call exit(1)
        endif

        ! Check that the maximum degree of the graph has been calculated
        ! properly by the initialization routine
        if (g%max_degree/=6) then
            print *, 'Max degree of g should be 6; it is:',g%max_degree
            call exit(1)
        endif

        ! Add in new edges to symmetrize the graph
        do i=1,12
            call g%add_edge(edges(2,i),edges(1,i))
        enddo

        ! Check that edge 1 is still connected to edge 2 and that edge 2
        ! is now connected to edge 1
        if (.not.g%connected(1,2) .or. .not.g%connected(2,1)) then
            print *, 'Supposed to have 1 -> 2 and 2 -> 1'
            call exit(1)
        endif

        ! Check that the total number of edges has been properly updated
        if (g%ne/=24) then
            print *, '12 edges were added; graph should have 24 edges'
            print *, 'Graph has ',g%ne,' edges'
            print *, test
            call exit(1)
        endif

        allocate(nbrs(g%max_degree))

        ! Check that finding all neighbors of a given edge works
        call g%neighbors(1,nbrs)
        correct = .true.
        do i=1,g%max_degree
            if (nbrs(i)/=i+1) correct = .false.
        enddo
        if (.not.correct) then
            print *, 'Node 1 should neighbor all other nodes'
            call exit(1)
        endif

        ! Check that iterating through a graph's edges works
        cursor = g%make_cursor(0)
        edg = g%get_edges(cursor,12,num_returned)

        found_by_iterator = .false.
        do k=1,12
            i = edg(1,k)
            j = edg(2,k)
            do l=1,12
                if (i==edges(1,l) .and. j==edges(2,l)) then
                    found_by_iterator(l) = .true.
                elseif (i==edges(2,l) .and. j==edges(1,l)) then
                    found_by_iterator(l+12) = .true.
                endif
            enddo
        enddo

        edg = g%get_edges(cursor,12,num_returned)

        do k=1,12
            i = edg(1,k)
            j = edg(2,k)
            do l=1,12
                if (i==edges(1,l) .and. j==edges(2,l)) then
                    found_by_iterator(l) = .true.
                elseif (i==edges(2,l) .and. j==edges(1,l)) then
                    found_by_iterator(l+12) = .true.
                endif
            enddo
        enddo

        correct = .true.
        do k=1,24
            correct = correct .and. found_by_iterator(k)
        enddo

        if (.not.correct) then
            print *, 'Iterating through graph edges failed'
            call exit(1)
        endif


        ! Check that permutation works properly
        allocate(p(7))
        do i=1,7
            p(i) = i-1
        enddo
        p(1) = 7
        call g%right_permute(p)
        call g%left_permute(p)
        deallocate(p)

        ! Delete all connections between node 1 and any even nodes
        do i=2,6,2
            call g%delete_edge(1,i)
            call g%delete_edge(i,1)
        enddo

        ! Check that edge deletion worked properly
        do i=2,6,2
            if (g%connected(1,i) .or. g%connected(i,1)) then
                print *, 'Test',test
                print *, 'All connections between 1 and even nodes were '
                print *, 'deleted, and yet 1 is still connected '
                print *, 'to node ',i
                call exit(1)
            endif
        enddo

        call g%free()

        deallocate(g,nbrs)

    enddo




end program graph_tests

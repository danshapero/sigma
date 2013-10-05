program graph_tests

use fempack

implicit none

    class(graph), allocatable :: g
    integer :: i,test
    integer, allocatable :: edges(:,:), nbrs(:), p(:)
    logical :: correct

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


    !----------------------------------------------------------------------!
    ! Loop through and test each graph type                                !
    !----------------------------------------------------------------------!
    do test=1,4

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

        call g%init(7,7,edges)

        if (g%ne/=12) then
            print *, 'Graph should have 12 edges'
            call exit(1)
        endif

        if (.not.g%connected(1,2) .or. g%connected(2,1)) then
            print *, 'Graph should have 1 -> 2, 2 -/> 1; we have:'
            if (.not.g%connected(1,2)) then
                print *, ' 1 -/> 2'
            else
                print *, ' 2 -> 1'
            endif
            call exit(1)
        endif

        if (g%max_degree/=6) then
            print *, 'Max degree of g should be 6; it is:',g%max_degree
            call exit(1)
        endif

        do i=1,12
            call g%add_edge(edges(2,i),edges(1,i))
        enddo

        if (.not.g%connected(1,2) .or. .not.g%connected(2,1)) then
            print *, 'Supposed to have 1 -> 2 and 2 -> 1'
        endif

        if (g%ne/=24) then
            print *, '12 edges were added; graph should have 24 edges'
            print *, 'Graph has ',g%ne,' edges'
            print *, test
        endif

        allocate(nbrs(g%max_degree))

        call g%neighbors(1,nbrs)
        correct = .true.
        do i=1,g%max_degree
            if (nbrs(i)/=i+1) correct = .false.
        enddo
        if (.not.correct) print *, 'Node 1 should neighbor all other nodes'

        allocate(p(7))
        do i=1,7
            p(i) = i-1
        enddo
        p(1) = 7
        call g%right_permute(p)
        call g%left_permute(p)
        deallocate(p)

        call g%free()

        deallocate(g,nbrs)

    enddo




end program graph_tests

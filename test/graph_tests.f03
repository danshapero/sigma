program graph_tests

use coo_graphs
use ll_graphs
use cs_graphs

implicit none

    class(graph), allocatable :: g
    integer :: i,j
    integer, allocatable :: edges(:,:)

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
    ! Test ll_graph type                                                   !
    !----------------------------------------------------------------------!
    allocate(ll_graph::g)

    call g%init(7,7,edges)

    if (.not.g%connected(1,2) .or. g%connected(2,1)) then
        print *, 'Graph should have 1 -> 2, 2 -/> 1; we have:'
        if (.not.g%connected(1,2)) then
            print *, ' 1 -/> 2'
        else
            print *, ' 2 -> 1'
        endif
    endif

    if (g%max_degree/=6) then
        print *, 'Max degree of g should be 6; value found:',g%max_degree
    endif

    do i=1,12
        call g%add_edge(edges(2,i),edges(1,i))
    enddo

    if (.not.g%connected(1,2) .or. .not.g%connected(2,1)) then
        print *, 'Supposed to have 1 -> 2 and 2 -> 1'
    endif



    call g%free()




    !----------------------------------------------------------------------!
    ! Test coo graph type                                                  !
    !----------------------------------------------------------------------!



    !----------------------------------------------------------------------!
    ! Test compressed sparse graph type                                    !
    !----------------------------------------------------------------------!



end program graph_tests

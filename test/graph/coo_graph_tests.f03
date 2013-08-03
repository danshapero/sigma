program coo_graph_tests

use coo_graphs

implicit none

    type(coo_graph) :: g
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

    call g%init(7,7,edges)

    print *, g%connected(1,2), g%connected(2,1)

    do i=1,12
        call g%add_edge(edges(2,i),edges(1,i))
    enddo

    print *, g%connected(1,2), g%connected(2,1)

    do i=2,7
        call g%delete_edge(1,i)
        call g%delete_edge(i,1)
    enddo

    print *, g%max_degree

end program coo_graph_tests

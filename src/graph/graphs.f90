module graphs

use graph_interface
use ll_graphs
use coo_graphs
use cs_graphs
use ellpack_graphs

implicit none


contains




!--------------------------------------------------------------------------!
subroutine graph_union(g,h1,h2,trans1,trans2)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    class(graph), intent(in) :: h1, h2
    logical, intent(in), optional :: trans1, trans2
    ! local variables
    integer :: i,j,k,n,order1(2),order2(2),nv1(2),nv2(2)
    integer :: num_blocks,num_returned,edges(2,64)
    type(graph_edge_cursor) :: cursor
    ! local versions of optional arguments
    logical :: tr1, tr2

    ! Get optional arguments
    tr1 = .false.
    tr2 = .false.
    if (present(trans1)) tr1 = trans1
    if (present(trans2)) tr2 = trans2

    ! Get the desired order for all the graphs
    order1 = [1, 2]
    order2 = [1, 2]
    if (tr1) order1 = [2, 1]
    if (tr2) order2 = [2, 1]

    ! Get the graph dimensions
    nv1 = [h1%n, h1%m]
    nv2 = [h2%n, h2%m]
    nv1 = nv1(order1)
    nv2 = nv2(order2)


    !------------------------------------------
    ! Check that the dimensions are consistent
    if (nv1(1)/=nv2(1) .or. nv1(2)/=nv2(2)) then
        print *, 'Inconsistent dimensions for graph union:'
        print *, 'graph 1:',nv1
        print *, 'graph 2:',nv2
        print *, 'Terminating.'
        call exit(1)
    endif


    !-------------------------------
    ! Initialize the output graph g
    call g%init(nv1(1),nv1(2),h1%max_degree+h2%max_degree)


    !---------------------------------------------
    ! Add the edges from h1 to the output graph g
    cursor = h1%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1
    do n=1,num_blocks
        edges = h1%get_edges(cursor,64,num_returned)

        do k=1,num_returned
            i = edges(order1(1),k)
            j = edges(order1(2),k)

            if (i/=0 .and. j/=0) call g%add_edge(i,j)
        enddo
    enddo


    !---------------------------------------------
    ! Add the edges from h2 to the output graph g
    cursor = h2%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1
    do n=1,num_blocks
        edges = h2%get_edges(cursor,64,num_returned)

        do k=1,num_returned
            i = edges(order2(1),k)
            j = edges(order2(2),k)

            if (i/=0 .and. j/=0) call g%add_edge(i,j)
        enddo
    enddo

end subroutine graph_union



end module graphs

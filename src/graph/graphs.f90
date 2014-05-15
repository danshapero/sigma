module graphs

use graph_interface
use ll_graphs
use coo_graphs
use cs_graphs
use ellpack_graphs

implicit none

private :: graph_product_optimized

contains




!--------------------------------------------------------------------------!
subroutine graph_union(g,h1,h2,trans1,trans2)                              !
!--------------------------------------------------------------------------!
!     Compute the union of two graphs h1, h2 and store the result in a new !
! graph g. This is used when computed the sum of two matrices.             !
!     The orientation of the edges of the input graphs h1, h2 could be     !
! reversed relative to the orientation of the output graph g, for example  !
! when adding a matrix in row-major ordering to a matrix in column-major   !
! ordering. In that case, an edge (v1,v2) of graph h1 is instead read as   !
! the edge (v2,v1).                                                        !
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
        call h1%get_edges(edges,cursor,64,num_returned)

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
        call h2%get_edges(edges,cursor,64,num_returned)

        do k=1,num_returned
            i = edges(order2(1),k)
            j = edges(order2(2),k)

            if (i/=0 .and. j/=0) call g%add_edge(i,j)
        enddo
    enddo

end subroutine graph_union



!--------------------------------------------------------------------------!
subroutine graph_product(g,h1,h2,trans1,trans2)                            !
!--------------------------------------------------------------------------!
!     Compute the product of two graphs h1, h2 and store the result in a   !
! third graph g. This is used when computing the product of two matrices.  !
!     The product of two graphs h1, h2 is a new graph such that (i,j) are  !
! connected in h1*h2 if there is a vertex k such that (i,k) are connected  !
! in h1 and (k,j) are connected in h2.                                     !
!     The orientation of the edges could be reversed, for the same reason  !
! reason that they could be reversed when computing a graph sum.           !
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
    ! an extra graph, for which querying all neighbors occurs in O(degree)
    type(ll_graph) :: hp


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
    if (nv1(2)/=nv2(1)) then
        print *, 'Inconsistent dimensions for graph product:'
        print *, 'graph 1:',nv1
        print *, 'graph 2:',nv2
        print *, 'Terminating.'
        call exit(1)
    endif


    ! Check and see if one of the graphs is in a nice enough format that
    ! the `neighbors` query takes O(degree) operations.

    ! If that's the case, we can use the `nice` version of the algorithm
    ! defined below.
    call graph_product_optimized(g,h1,h2, &
        & trans_h1=tr1, trans_h2=tr2, trans_g=.false.)

    ! If it's not the case, we will, at the expense of extra memory usage,
    ! make a copy of one of the graphs which is in a convenient format.

    ! Then go to the nice version of the algorithm.

end subroutine graph_product



!--------------------------------------------------------------------------!
subroutine graph_product_optimized(g,h1,h2,trans_h1,trans_h2,trans_g)      !
!--------------------------------------------------------------------------!
!     Compute the product of the graphs h1, h2 and store them in a third   !
! graph g.                                                                 !
!     This is an optimized version of the algorithm which assumes that the !
! graph h1 is in a format where the `neighbors` operation runs fast. The   !
! "user-facing" algorithm graph_product defers to this private, optimized  !
! version, either in the case that one of the graphs is already in a nice  !
! format, or in case it copies one of the graphs to a nice format and then !
! calls the optimized version.                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    class(graph), intent(in) :: h1, h2
    logical, intent(in), optional :: trans_h1, trans_h2, trans_g
    ! local variables
    type(ll_graph) :: gl
    integer :: i,j,k,l,m,n,d,order1(2),order2(2),nv1(2),nv2(2),ind(2)
    integer :: num_blocks,num_returned,edges(2,64)
    type(graph_edge_cursor) :: cursor
    ! a neighbors array
    integer, allocatable :: neighbors(:)
    ! local versions of optional arguments
    logical :: tr1, tr2, trg

    ! Get optional arguments
    tr1 = .false.
    tr2 = .false.
    trg = .false.
    if (present(trans_h1)) tr1 = trans_h1
    if (present(trans_h2)) tr2 = trans_h2
    if (present(trans_g)) trg = trans_g

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

    allocate(neighbors(h2%max_degree))


    !-----------------------------------------------
    ! Make a temporary graph for storing the result
    call gl%init(nv1(1),nv2(2))


    !-------------------------------------
    ! Iterate through all the edges of h1
    cursor = h1%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    do n=1,num_blocks
        call h1%get_edges(edges,cursor,64,num_returned)
        do l=1,num_returned
            ind = edges(order1,l)

            if (ind(1)/=0 .and. ind(2)/=0) then
                i = ind(order2(1))
                k = ind(order2(2))

                call h2%get_neighbors(neighbors,k)
                d = h2%degree(k)
                do m=1,d
                    j = neighbors(m)
                    call gl%add_edge(i,j)
                enddo
            endif
        enddo
    enddo



    !------------------------------------------------
    ! Copy the linked list graph to the output graph
    call g%copy(gl,trg)
    call gl%free()

end subroutine graph_product_optimized



end module graphs

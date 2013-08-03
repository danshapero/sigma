module cs_graphs

use abstract_graphs
use util

implicit none



!--------------------------------------------------------------------------!
type, extends(graph) :: cs_graph                                           !
!--------------------------------------------------------------------------!
    integer, allocatable :: ia(:), ja(:)
contains
    procedure :: init => cs_init
    procedure :: neighbors => cs_neighbors
    procedure :: connected => cs_connected
    procedure :: find_edge => cs_find_edge
    procedure :: add_edge  => cs_add_edge
    procedure :: delete_edge => cs_delete_edge
    procedure :: free => cs_free
    procedure :: dump_edges => cs_dump_edges
    ! auxiliary routines
    procedure :: sort_ja
end type cs_graph





contains



!--------------------------------------------------------------------------!
subroutine cs_init(g,n,m,edges)                                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, edges(:,:)
    ! local variables
    integer :: i,k,work(n),ne

    ! Set the number of (left-)nodes in the graph and allocate the node
    ! pointer array ia
    g%n = n
    allocate( g%ia(n+1) )

    ! If the graph is not square, set the number of right-nodes
    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    ! If the graph is being initialized with a set of edges, set the number
    ! of edges and allocate space in the node array ja
    if (present(edges)) then
        ne = size(edges,2)
        g%ne = ne
    else
        ne = 0
        g%ne = ne
    endif
    allocate( g%ja(ne) )

    if (present(edges)) then
        work = 0
        do k=1,ne
            i = edges(1,k)
            work(i) = work(i)+1
        enddo

        g%max_degree = maxval(work)

        g%ia(1) = 1
        do k=2,n+1
            g%ia(k) = g%ia(k-1)+work(k-1)
        enddo

        work = 0
        do k=1,ne
            i = edges(1,k)
            g%ja( g%ia(i)+work(i) ) = edges(2,k)
            work(i) = work(i)+1
        enddo

        call g%sort_ja()
    else
        g%max_degree = 0
    endif

end subroutine cs_init



!--------------------------------------------------------------------------!
subroutine cs_neighbors(g,i,nbrs)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(in) :: i
    integer, intent(out) :: nbrs(:)
    ! local variables
    integer :: start,finish,degree

    degree = g%ia(i+1)-g%ia(i)
    start = g%ia(i)
    finish = g%ia(i+1)-1

    nbrs = 0
    nbrs(1:degree) = g%ja(start:finish)

end subroutine cs_neighbors



!--------------------------------------------------------------------------!
function cs_connected(g,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: cs_connected
    ! local variables
    integer :: k

    cs_connected = .false.

    do k=g%ia(i),g%ia(i+1)-1
        if (g%ja(k)==j) cs_connected = .true.
    enddo

end function cs_connected



!--------------------------------------------------------------------------!
function cs_find_edge(g,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(in) :: i,j
    integer :: cs_find_edge
    ! local variables
    integer :: k

    cs_find_edge = g%ia(i)

    do k=g%ia(i),g%ia(i+1)-1
        if (g%ja(k)<j) then
            cs_find_edge = k+1
        endif
    enddo

end function cs_find_edge



!--------------------------------------------------------------------------!
subroutine cs_add_edge(g,i,j)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k,indx,ja_temp(g%ne)

    if (.not.g%connected(i,j)) then
        ! Find the index in the list of edges where the new edge is to be
        ! inserted
        indx = g%find_edge(i,j)

        ! Copy the current structure into a temporary one
        ja_temp = g%ja

        ! Deallocate the old structure and make a new one with extra room
        deallocate(g%ja)
        allocate(g%ja(g%ne+1))

        ! Build the new array ja
        g%ja(1:indx-1) = ja_temp(1:indx-1)
        g%ja(indx) = j
        g%ja(indx+1:g%ne+1) = ja_temp(indx:g%ne)

        ! The number of edges in the graph is incremented by 1, and the
        ! starting index in ja is incremented for all nodes after the node
        ! into which the new edge is to be inserted
        g%ne = g%ne+1
        do k=i+1,g%n+1
            g%ia(k) = g%ia(k)+1
        enddo

        ! If the degree of node i is now the greatest of all nodes in the
        ! graph, update the degree accordingly
        if (g%ia(i+1)-g%ia(i)>g%max_degree) then
            g%max_degree = g%ia(i+1)-g%ia(i)
        endif
    endif

end subroutine cs_add_edge



!--------------------------------------------------------------------------!
subroutine cs_delete_edge(g,i,j)                                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k,indx,degree,ja_temp(g%ne)

    if (g%connected(i,j)) then
        ! Record the degree of the node from which an edge is to be
        ! removed; we use this later to check whether the maximum degree
        ! of the graph has decreased
        degree = g%ia(i+1)-g%ia(i)

        ! Find the index in the list of edges of the edge to be removed
        indx = g%find_edge(i,j)

        ja_temp = g%ja

        deallocate(g%ja)
        allocate(g%ja(g%ne-1))

        ! Build the new array ja
        g%ja(1:indx-1) = ja_temp(1:indx-1)
        g%ja(indx:g%ne-1) = ja_temp(indx+1:g%ne)

        ! The number of edges in the graph is decremented by 1, and the
        ! starting index in ja is decremented for all nodes after the node
        ! from which the edge was removed
        g%ne = g%ne-1
        do k=i+1,g%n
            g%ia(k) = g%ia(k)-1
        enddo

        ! If the degree of node i was the greatest among all nodes in the
        ! entire graph, check to see if max_degree needs to be decremented
        if (degree==g%max_degree) then
            ! Tentatively set the maximum degree of the graph to the (now 
            ! smaller) degree of the node i from which the edge was removed.
            g%max_degree = g%ia(i+1)-g%ia(i)

            ! Iterate through all the nodes;
            do k=1,g%n
                ! if any node has a degree greater than the tentative
                ! max_degree, then increment the maximum degree.
                degree = g%ia(k+1)-g%ia(k)
                if (degree>g%max_degree) g%max_degree = degree
            enddo
        endif
    endif

end subroutine cs_delete_edge



!--------------------------------------------------------------------------!
subroutine cs_free(g)                                                      !
!--------------------------------------------------------------------------!
    class(cs_graph), intent(inout) :: g

    deallocate(g%ia,g%ja)
    g%n = 0
    g%m = 0
    g%ne = 0
    g%max_degree = 0

end subroutine cs_free



!--------------------------------------------------------------------------!
subroutine cs_dump_edges(g,edges)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: i,k

    do i=1,g%n
        do k=g%ia(i),g%ia(i+1)-1
            edges(1,k) = i
            edges(2,k) = g%ja(k)
        enddo
    enddo

end subroutine cs_dump_edges



!--------------------------------------------------------------------------!
subroutine sort_ja(g)                                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    ! local variables
    integer :: k,start,finish,num,p(g%max_degree)

    do k=1,g%n
        start  = g%ia(k)
        finish = g%ia(k+1)-1
        num = finish-start+1
        p(1:num) = order(g%ja(start:finish))
        g%ja(start:finish) = g%ja( p(1:num) )
    enddo

end subroutine sort_ja



end module cs_graphs

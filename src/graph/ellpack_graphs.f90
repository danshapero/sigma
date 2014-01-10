module ellpack_graphs

use graphs
use util

implicit none


!--------------------------------------------------------------------------!
type, extends(graph) :: ellpack_graph                                      !
!--------------------------------------------------------------------------!
    integer, allocatable :: node(:,:)
contains
    procedure :: init => ellpack_graph_init
    procedure :: neighbors => ellpack_neighbors
    procedure :: connected => ellpack_connected
    procedure :: find_edge => ellpack_find_edge
    procedure :: make_cursor => ellpack_make_cursor
    procedure :: get_edges => ellpack_get_edges
    procedure :: add_edge => ellpack_add_edge
    procedure :: delete_edge => ellpack_delete_edge
    procedure :: left_permute => ellpack_graph_left_permute
    procedure :: right_permute => ellpack_graph_right_permute
    procedure :: free => ellpack_free
    procedure :: dump_edges => ellpack_dump_edges
    ! auxiliary routines
    procedure :: ellpack_add_edge_with_max_degree_increase
    procedure :: ellpack_max_degree_decrease
end type ellpack_graph





contains



!--------------------------------------------------------------------------!
subroutine ellpack_graph_init(g,n,m,edges)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, edges(:,:)
    ! local variables
    integer :: i,j,k,ne,degree(n)

    g%n = n

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    if (present(edges)) then
        ne = size(edges,2)
        g%ne = ne

        degree = 0
        do k=1,ne
            i = edges(1,k)
            degree(i) = degree(i)+1
        enddo
        g%max_degree = maxval(degree)
    else
        ne = 0
        g%ne = ne
        g%max_degree = 0
    endif

    allocate( g%node(g%max_degree,g%n) )
    g%node = 0

    if (present(edges)) then
        degree = 0

        do k=1,ne
            i = edges(1,k)
            j = edges(2,k)

            degree(i) = degree(i)+1
            g%node(degree(i),i) = j
        enddo
    endif

end subroutine ellpack_graph_init



!--------------------------------------------------------------------------!
subroutine ellpack_neighbors(g,i,nbrs)                                     !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i
    integer, intent(out) :: nbrs(:)

    nbrs(1:g%max_degree) = g%node(:,i)

end subroutine ellpack_neighbors



!--------------------------------------------------------------------------!
function ellpack_connected(g,i,j)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: ellpack_connected
    ! local variables
    integer :: k

    ellpack_connected = .false.

    do k=1,g%max_degree
        if (g%node(k,i)==j) ellpack_connected = .true.
    enddo

end function ellpack_connected



!--------------------------------------------------------------------------!
function ellpack_find_edge(g,i,j)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i,j
    integer :: ellpack_find_edge
    ! local variables
    integer :: k,l,total,degree

    ellpack_find_edge = -1

    total = 0
    do l=1,i-1
        degree = g%max_degree
        do k=g%max_degree,1,-1
            if (g%node(k,l)==0) then
                degree = k-1
            endif
        enddo

        total = total+degree
    enddo

    do k=1,g%max_degree
        if (g%node(k,i)==j) ellpack_find_edge = total+k
    enddo

    !! This approach is quite blunt and we need to do something else if we
    !! we want to keep ellpack graphs as a performant rather than an easy
    !! format.

end function ellpack_find_edge



!--------------------------------------------------------------------------!
function ellpack_make_cursor(g,thread) result(cursor)                      !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: thread
    type(graph_edge_cursor) :: cursor

    cursor%start = 1
    cursor%final = g%ne
    cursor%current = 0
    cursor%edge = [1, g%node(1,1)]
    cursor%indx = 0

end function ellpack_make_cursor



!--------------------------------------------------------------------------!
function ellpack_get_edges(g,cursor,num_edges,num_returned) result(edges)  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    integer :: edges(2,num_edges)
    ! local variables
    integer :: i,k,degree,num_added,num_from_this_row

    ! Set up the returned edges to be 0
    edges = 0

    ! Find out how many nodes' edges the current request encompasses
    num_returned = min(num_edges,cursor%final-cursor%current)
    i = cursor%edge(1)

    ! Loop over each of those nodes
    num_added = 0
    do while(num_added<num_returned)
        ! Compute the degree of node i
        degree = g%max_degree
        do k=g%max_degree,1,-1
            if (g%node(k,i)==0) then
                degree = k-1
            endif
        enddo

        ! Find how many nodes to return from the current row
        num_from_this_row = min(degree-cursor%indx, num_returned-num_added)

        do k=1,num_from_this_row
            edges(1,num_added+k) = i
            edges(2,num_added+k) = g%node(cursor%indx+k,i)
        enddo

        cursor%indx = mod(cursor%indx+num_from_this_row,degree)

        ! If we returned all nodes from this row, increment the row
        if (num_from_this_row == degree-cursor%indx) then
            i = i+1
        endif

        ! Increase the number of edges added
        num_added = num_added+num_from_this_row
    enddo

    cursor%current = cursor%current+num_returned
    cursor%edge(1) = i

end function ellpack_get_edges



!--------------------------------------------------------------------------!
subroutine ellpack_add_edge(g,i,j)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k,indx

    ! If the two nodes i,j are already connected, we needn't add the edge
    if (.not.g%connected(i,j)) then
        ! Find the index in g%node(:,i) where we can add in node j
        indx = -1

        do k=1,g%max_degree
            if (g%node(k,i)==0) then
                indx = k
                exit
            endif
        enddo

        ! If there is room to add j, then do so
        if (indx/=-1) then
            g%node(indx,i) = j

        ! If there is no room, that means that degree(i) = max degree of g.
        ! We then have to rebuild the array g%node to reflect the fact that
        ! the graph has a higher max degree.
        else
            call g%ellpack_add_edge_with_max_degree_increase(i,j)
        endif

        g%ne = g%ne+1
    endif

end subroutine ellpack_add_edge



!--------------------------------------------------------------------------!
subroutine ellpack_delete_edge(g,i,j)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k,indx
    logical :: decrease_max_degree

    ! If nodes i,j are not connected to begin with, there is no edge to
    ! delete and thus nothing to do
    if (g%connected(i,j)) then
        indx = g%find_edge(i,j)

        g%node(indx:g%max_degree-1,i) = g%node(indx+1:g%max_degree,i)
        g%node(g%max_degree,i) = 0

        decrease_max_degree = .true.
        do k=1,g%n
            if (g%node(g%max_degree,k)/=0) then
                decrease_max_degree = .false.
                exit
            endif
        enddo

        if (decrease_max_degree) call g%ellpack_max_degree_decrease()

        g%ne = g%ne-1
    endif

end subroutine ellpack_delete_edge



!--------------------------------------------------------------------------!
subroutine ellpack_graph_left_permute(g,p)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,node(g%max_degree,g%n)

    do i=1,g%n
        node(:,p(i)) = g%node(:,i)
    enddo

    g%node = node

end subroutine ellpack_graph_left_permute



!--------------------------------------------------------------------------!
subroutine ellpack_graph_right_permute(g,p)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,j,k

    do i=1,g%n
        do k=1,g%max_degree
            j = g%node(k,i)
            if (j/=0) g%node(k,i) = p(j)
        enddo
    enddo

end subroutine ellpack_graph_right_permute



!--------------------------------------------------------------------------!
subroutine ellpack_free(g)                                                 !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(inout) :: g

    deallocate(g%node)
    g%n = 0
    g%m = 0
    g%ne = 0
    g%max_degree = 0

end subroutine ellpack_free



!--------------------------------------------------------------------------!
subroutine ellpack_dump_edges(g,edges)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: i,j,k,next

    next = 0
    do i=1,g%n
        do k=1,g%max_degree
            j = g%node(k,i)
            if (j/=0) then
                next = next+1
                edges(1,next) = i
                edges(2,next) = j
            else
                exit
            endif
        enddo
    enddo

end subroutine ellpack_dump_edges



!--------------------------------------------------------------------------!
subroutine ellpack_add_edge_with_max_degree_increase(g,i,j)                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k, node(g%max_degree,g%n)

    node = g%node
    deallocate(g%node)
    allocate(g%node(g%max_degree+1,g%n))

    g%node = 0
    do k=1,g%n
        g%node(1:g%max_degree,k) = node(1:g%max_degree,k)
    enddo

    g%node(g%max_degree+1,i) = j

    g%max_degree = g%max_degree+1

end subroutine ellpack_add_edge_with_max_degree_increase



!--------------------------------------------------------------------------!
subroutine ellpack_max_degree_decrease(g)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    ! local variables
    integer :: k, node(g%max_degree,g%n)

    node = g%node
    deallocate(g%node)
    allocate(g%node(g%max_degree-1,g%n))

    g%node = 0
    do k=1,g%n
        g%node(1:g%max_degree-1,k) = node(1:g%max_degree-1,k)
    enddo

    g%max_degree = g%max_degree-1

end subroutine ellpack_max_degree_decrease





end module ellpack_graphs

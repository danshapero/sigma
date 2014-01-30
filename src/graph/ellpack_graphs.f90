module ellpack_graphs

use graphs
use util

implicit none


!--------------------------------------------------------------------------!
type, extends(graph) :: ellpack_graph                                      !
!--------------------------------------------------------------------------!
    integer, allocatable :: node(:,:)
    integer :: max_neighbors
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
    procedure :: left_permute_edge_reorder => ellpack_left_perm_edge_reorder
    procedure :: right_permute_edge_reorder &
                                        & => ellpack_right_perm_edge_reorder
    procedure :: free => ellpack_free
    procedure :: dump_edges => ellpack_dump_edges
    ! auxiliary routines
    procedure :: ellpack_add_edge_with_max_degree_increase
    procedure :: ellpack_max_degree_decrease
end type ellpack_graph





contains



!--------------------------------------------------------------------------!
subroutine ellpack_graph_init(g,n,m,num_neighbor_nodes)                    !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, num_neighbor_nodes(:)

    g%n = n

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    if (present(num_neighbor_nodes)) then
        g%max_neighbors = maxval(num_neighbor_nodes)
    else
        g%max_neighbors = 6
    endif
    allocate(g%node(g%max_neighbors,g%n))
    g%node = 0
    g%max_degree = 0
    g%ne = 0
    g%capacity = g%max_neighbors*g%n

end subroutine ellpack_graph_init



!--------------------------------------------------------------------------!
subroutine ellpack_neighbors(g,i,nbrs)                                     !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i
    integer, intent(out) :: nbrs(:)

    nbrs = 0
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
    integer :: k

    ellpack_find_edge = -1

    do k=g%max_degree,1,-1
        if (g%node(k,i)==j) ellpack_find_edge = (i-1)*g%max_degree+k
    enddo

end function ellpack_find_edge



!--------------------------------------------------------------------------!
function ellpack_make_cursor(g,thread) result(cursor)                      !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: thread
    type(graph_edge_cursor) :: cursor

    cursor%start = 1
    cursor%final = g%capacity
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
    integer :: i,i1,i2,num_added,num_from_this_row

    ! Set up the returned edges to be 0
    edges = 0

    ! Find out how many nodes' edges the current request encompasses
    num_returned = min(num_edges,cursor%final-cursor%current)

    ! Find the starting and ending nodes for this edge retrieval
    i1 = cursor%current/g%max_neighbors+1
    i2 = (cursor%current+num_returned-1)/g%max_neighbors+1

    ! Set the number of edges added to 0
    num_added = 0

    ! Loop from the starting node to the ending node
    do i=i1,i2
        ! Find how many edges we're retrieving from this row
        num_from_this_row = min(g%max_degree-cursor%indx, &
                                & num_returned-num_added)

        ! Fill in the return array
        edges(1,num_added+1:num_added+num_from_this_row) = i
        edges(2,num_added+1:num_added+num_from_this_row) = &
            & g%node(cursor%indx+1:cursor%indx+num_from_this_row,i)

        ! Increment the number of edges added
        num_added = num_added+num_from_this_row

        ! Modify the index storing the place within the row that we left off
        cursor%indx = mod(cursor%indx+num_from_this_row,g%max_degree)
    enddo

    cursor%current = cursor%current+num_returned

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

        do k=g%max_neighbors,1,-1
            if (g%node(k,i)==0) then
                indx = k
            endif
        enddo

        ! If there is room to add j, then do so
        if (indx/=-1) then
            g%node(indx,i) = j
            g%ne = g%ne+1
        ! If there is no room, that means that degree(i) = max degree of g.
        else
            print *, 'Not enough space to add edge',i,j
        endif

        ! Increment the maximum degree of the graph if need be
        g%max_degree = max(g%max_degree,indx)
    endif

end subroutine ellpack_add_edge



!--------------------------------------------------------------------------!
subroutine ellpack_delete_edge(g,i,j)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: indx
    logical :: max_degree_decrease

    ! If nodes i,j are not connected to begin with, there is no edge to
    ! delete and thus nothing to do
    if (g%connected(i,j)) then
        ! Set a boolean to be true if we're removing an edge from a node
        ! of maximum degree
        max_degree_decrease = (g%node(g%max_degree,i)/=0)

        ! Find the location indx in memory where edge (i,j) is stored
        do indx=1,g%max_degree
            if (g%node(indx,i)==j) exit
        enddo

        ! Overwrite indx with the other nodes connected to i
        g%node(indx:g%max_degree-1,i) = g%node(indx+1:g%max_degree,i)

        ! Zero out the last node connected to i
        g%node(g%max_degree,i) = 0

        ! Decrement the number of edges
        g%ne = g%ne-1

        ! If node i had maximum degree, check all the other nodes to see if
        ! the max degree has decreased
        if (max_degree_decrease) then
            max_degree_decrease = (maxval(g%node(g%max_degree,:))==0)
        endif

        ! If so, decrement the max_degree member of g
        if (max_degree_decrease) g%max_degree = g%max_degree-1
    endif

end subroutine ellpack_delete_edge



!--------------------------------------------------------------------------!
subroutine ellpack_graph_left_permute(g,p)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,node(g%max_neighbors,g%n)

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
subroutine ellpack_left_perm_edge_reorder(g,p,edge_p)                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out) :: edge_p(:,:)
    ! local variables
    integer :: i

    ! Permute the graph
    call g%left_permute(p)

    ! Report the resulting edge permutation
    allocate(edge_p(3,g%n))

    do i=1,g%n
        edge_p(1,i) = g%max_neighbors*(i-1)+1
        edge_p(2,i) = g%max_neighbors*(p(i)-1)+1
        edge_p(3,i) = g%max_neighbors
    enddo

end subroutine ellpack_left_perm_edge_reorder



!--------------------------------------------------------------------------!
subroutine ellpack_right_perm_edge_reorder(g,p,edge_p)                     !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out) :: edge_p(:,:)

    call g%right_permute(p)
    allocate(edge_p(0,0))

end subroutine ellpack_right_perm_edge_reorder



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

module cs_graphs

use graphs
use util

implicit none



!--------------------------------------------------------------------------!
type, extends(graph) :: cs_graph                                           !
!--------------------------------------------------------------------------!
    integer, allocatable :: ptr(:), node(:), degree(:)
contains
    procedure :: init => cs_init
    procedure :: neighbors => cs_neighbors
    procedure :: connected => cs_connected
    procedure :: find_edge => cs_find_edge
    procedure :: make_cursor => cs_make_cursor
    procedure :: get_edges => cs_get_edges
    procedure :: add_edge  => cs_add_edge
    procedure :: delete_edge => cs_delete_edge
    procedure :: left_permute => cs_graph_left_permute
    procedure :: right_permute => cs_graph_right_permute
    procedure :: free => cs_free
    procedure :: dump_edges => cs_dump_edges
    ! auxiliary routines
    procedure :: sort_node
    procedure, private :: max_degree_update
end type cs_graph





contains



!--------------------------------------------------------------------------!
subroutine cs_init(g,n,m,num_neighbor_nodes)                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, num_neighbor_nodes(:)
    ! local variables
    integer :: k,ne

    ! Set the number of (left-)nodes in the graph and allocate the node
    ! pointer array ptr
    g%n = n
    allocate( g%ptr(n+1) )

    ! If the graph is not square, set the number of right-nodes
    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    ! If the graph is being initialized with a set of edges, set the number
    ! of edges and allocate space in the node array node
    if (present(num_neighbor_nodes)) then
        ne = sum(num_neighbor_nodes)
    else
        ne = g%n
    endif

    allocate( g%node(ne) )
    g%node = 0
    g%capacity = ne

    if (present(num_neighbor_nodes)) then
        g%ptr(1) = 1
        do k=1,g%n
            g%ptr(k+1) = g%ptr(k)+num_neighbor_nodes(k)
        enddo
    else
        g%ptr = 1
    endif

    g%ne = 0
    g%max_degree = 0

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

    degree = g%ptr(i+1)-g%ptr(i)
    start = g%ptr(i)
    finish = g%ptr(i+1)-1

    nbrs = 0
    nbrs(1:degree) = g%node(start:finish)

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

    do k=g%ptr(i),g%ptr(i+1)-1
        if (g%node(k)==j) cs_connected = .true.
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

    cs_find_edge = -1

    do k=g%ptr(i),g%ptr(i+1)-1
        if (g%node(k)==j) cs_find_edge = k
    enddo

end function cs_find_edge



!--------------------------------------------------------------------------!
function cs_make_cursor(g,thread) result(cursor)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(in) :: thread
    type(graph_edge_cursor) :: cursor
    ! local variables
    integer :: k

    cursor%start = 1
    cursor%final = g%capacity
    cursor%current = 0

    k = 1
    do while(g%ptr(k+1)-g%ptr(k)==0)
        k = k+1
    enddo

    cursor%edge = [k, g%node(1)]

end function cs_make_cursor



!--------------------------------------------------------------------------!
function cs_get_edges(g,cursor,num_edges,num_returned) result(edges)       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    integer :: edges(2,num_edges)
    ! local variables
    integer :: i, k

    ! Set up the returned edges to be 0
    edges = 0

    ! Count how many edges we're actually going to return; we'll either
    ! return how many edges the user asked for, or, if that amount would
    ! go beyond the final edge that the cursor is allowed to access, all
    ! of the remaining edges
    num_returned = min(num_edges,cursor%final-cursor%current)

    ! Fill the edges array's second row with the edge endpoints
    edges(2,1:num_returned) = &
        & g%node(cursor%current+1:cursor%current+num_returned)

    ! Fill in the edges array's first row with the edge start points
    i = cursor%edge(1)

    do k=1,num_returned
        !! This is going to produce code that cannot be analyzed and
        !! optimized well by the compiler. Would be good to find something
        !! more slick with a predictable access pattern.
        do while(cursor%current >= g%ptr(i+1)-1)
            i = i+1
        enddo

        edges(1,k) = i
        cursor%current = cursor%current+1
    enddo

    cursor%edge = edges(:,num_returned)

end function cs_get_edges



!--------------------------------------------------------------------------!
subroutine cs_add_edge(g,i,j)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k
    logical :: added

    if (.not.g%connected(i,j)) then
        ! Try to see if the new neighbor can be added without reallocating
        ! memory
        added = .false.
        do k=g%ptr(i),g%ptr(i+1)-1
            if (g%node(k)==0) then
                g%node(k) = j
                added = .true.
                exit
            endif
        enddo

        g%ne = g%ne+1

        if (.not.added) then
            print *, 'Not enough space to add edge',i,j
        endif

        if (k-g%ptr(i)+1>g%max_degree) then
            g%max_degree = k-g%ptr(i)+1
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
    integer :: jt,k,indx,degree

    ! Find the index in the list of edges of the edge to be removed
    indx = g%find_edge(i,j)

    if (indx/=-1) then
        ! Record the degree of the node from which an edge is to be
        ! removed, and find the last node that i is connected to
        do k=g%ptr(i),g%ptr(i+1)-1
            if (g%node(k)/=0) then
                degree = k-g%ptr(i)+1
                jt = g%node(k)
            endif
        enddo

        ! If there were more than two edges connected to node i, then
        ! replace the edge (i,j) with the removed edge (i,jt)
        if (degree>1) then
            g%node(indx) = jt
        endif

        ! Remove the last edge (i,jt) connected to i
        g%node( g%ptr(i)+degree-1 ) = 0

        if (degree==g%max_degree) then
            call g%max_degree_update()
        endif

        ! Decrement the number of edges in g
        g%ne = g%ne-1
    endif

end subroutine cs_delete_edge



!--------------------------------------------------------------------------!
subroutine cs_graph_left_permute(g,p)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,k,ptr(g%n+1),node(g%capacity)

    do i=1,g%n
        ptr(p(i)+1) = g%ptr(i+1)-g%ptr(i)
    enddo

    ptr(1) = 1
    do i=1,g%n
        ptr(i+1) = ptr(i+1)+ptr(i)
    enddo

    do i=1,g%n
        do k=0,g%ptr(i+1)-g%ptr(i)-1
            node( ptr(p(i))+k ) = g%node( g%ptr(i)+k )
        enddo
    enddo

    g%ptr = ptr
    g%node = node

end subroutine cs_graph_left_permute



!--------------------------------------------------------------------------!
subroutine cs_graph_right_permute(g,p)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    ! local variables
    integer :: j,k

    do k=1,g%capacity
        j = g%node(k)
        if (j/=0) g%node(k) = p(j)
    enddo

end subroutine cs_graph_right_permute



!--------------------------------------------------------------------------!
subroutine cs_free(g)                                                      !
!--------------------------------------------------------------------------!
    class(cs_graph), intent(inout) :: g

    deallocate(g%ptr,g%node)
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
        do k=g%ptr(i),g%ptr(i+1)-1
            edges(1,k) = i
            edges(2,k) = g%node(k)
        enddo
    enddo

end subroutine cs_dump_edges



!--------------------------------------------------------------------------!
subroutine sort_node(g)                                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    ! local variables
    integer :: k,start,finish,num,p(g%max_degree),node(g%max_degree)

    do k=1,g%n
        start  = g%ptr(k)
        finish = g%ptr(k+1)-1
        num = finish-start+1
        node(1:num) = g%node(start:finish)
        p(1:num) = order(node(1:num))
        g%node(start:finish) = node(p(1:num))
    enddo

end subroutine sort_node



!--------------------------------------------------------------------------!
subroutine max_degree_update(g)                                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    ! local variables
    integer :: i,k,degree,max_d

    max_d = 0

    ! loop through all the nodes
    do i=1,g%n
        ! set the degree of i to be 0
        degree = 0

        ! check how many neighbors i has
        do k=g%ptr(i),g%ptr(i+1)-1
            if (g%node(k)/=0) degree = k-g%ptr(i)+1
        enddo

        ! update the max degree of g if node i has a greater degree
        max_d = max(max_d,degree)
    enddo

    g%max_degree = max_d

end subroutine max_degree_update



!--------------------------------------------------------------------------!
subroutine cs_compress(g)                                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    ! local variables
    integer :: i,j,k,ptr(g%n+1),node(g%capacity)

    ptr = g%ptr
    node = g%node

    g%ptr = 0
    g%ptr(1) = 1
    g%ne = 0

    do i=1,g%n
        do k=ptr(i),ptr(i+1)-1
            j = node(k)
            if (j/=0) then
                g%ptr(i+1) = g%ptr(i+1)+1
                g%ne = g%ne+1
                g%node(g%ne) = j
            endif
        enddo

        g%ptr(i+1) = g%ptr(i)+g%ptr(i+1)
    enddo

end subroutine cs_compress


end module cs_graphs

module ll_graphs

use graphs
use types, only: dynamic_array

implicit none


!--------------------------------------------------------------------------!
type, extends(graph) :: ll_graph                                           !
!--------------------------------------------------------------------------!
    type(dynamic_array), allocatable :: lists(:)
contains
    procedure :: init => ll_init
    procedure :: neighbors => ll_neighbors
    procedure :: connected => ll_connected
    procedure :: find_edge => ll_find_edge
    procedure :: make_cursor => ll_make_cursor
    procedure :: get_edges => ll_get_edges
    procedure :: add_edge => ll_add_edge
    procedure :: delete_edge => ll_delete_edge
    procedure :: left_permute => ll_graph_left_permute
    procedure :: right_permute => ll_graph_right_permute
    procedure :: left_permute_edge_reorder => ll_left_perm_edge_reorder
    procedure :: right_permute_edge_reorder => ll_right_perm_edge_reorder
    procedure :: free => ll_free
    procedure :: dump_edges => ll_dump_edges
end type ll_graph





contains



!--------------------------------------------------------------------------!
subroutine ll_init(g,n,m,num_neighbor_nodes)                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, num_neighbor_nodes(:)
    ! local variables
    integer :: k,ne

    g%n = n
    allocate(g%lists(n))

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    if (present(num_neighbor_nodes)) then
        do k=1,n
            call g%lists(k)%init(capacity=num_neighbor_nodes(k), &
                & min_capacity=2)
        enddo
        ne = sum(num_neighbor_nodes)
    else
        do k=1,n
            call g%lists(k)%init(capacity=4, min_capacity=2)
        enddo
        ne = 4*g%n
    endif

    g%ne = 0
    g%max_degree = 0
    g%capacity = ne

end subroutine ll_init



!--------------------------------------------------------------------------!
subroutine ll_neighbors(g,i,nbrs)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: i
    integer, intent(out) :: nbrs(:)
    ! local variables
    integer :: k

    nbrs = 0
    do k=1,g%lists(i)%length
        nbrs(k) = g%lists(i)%get_entry(k)
    enddo

end subroutine ll_neighbors



!--------------------------------------------------------------------------!
function ll_connected(g,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: ll_connected
    ! local variables
    integer :: k

    ll_connected = .false.
    do k=1,g%lists(i)%length
        if (g%lists(i)%get_entry(k)==j) then
            ll_connected = .true.
            exit
        endif
    enddo

end function ll_connected



!--------------------------------------------------------------------------!
function ll_find_edge(g,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: i,j
    integer :: ll_find_edge
    ! local variables
    integer :: k,total

    ll_find_edge = -1

    total = 0
    do k=1,i-1
        total = total+g%lists(k)%length
    enddo

    do k=1,g%lists(i)%length
        if (g%lists(i)%get_entry(k)==j) then
            ll_find_edge = total+k
            exit
        endif
    enddo

end function ll_find_edge



!--------------------------------------------------------------------------!
function ll_make_cursor(g,thread) result(cursor)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: thread
    type(graph_edge_cursor) :: cursor
    ! local variables
    integer :: k

    cursor%start = 1
    cursor%final = g%ne
    cursor%current = 0

    k = 1
    do while (g%lists(k)%length==0)
        k = k+1
    enddo

    cursor%edge = [k, g%lists(k)%get_entry(1)]
    cursor%indx = 0

end function ll_make_cursor



!--------------------------------------------------------------------------!
function ll_get_edges(g,cursor,num_edges,num_returned) result(edges)       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    integer :: edges(2,num_edges)
    ! local variables
    integer :: i,k,num_added,num_from_this_row

    ! Set up the returned edges to be 0
    edges = 0

    ! Count how many edges we're actually going to return
    num_returned = min(num_edges, cursor%final-cursor%current)

    ! Find the last row we left off at
    i = cursor%edge(1)

    num_added = 0
    do while(num_added<num_returned)
        ! Either we are returning all the edges for the current node i,
        ! or the current request doesn't call for so many edges and we are
        ! returning fewer
        num_from_this_row = min(g%lists(i)%length-cursor%indx, &
                                & num_returned-num_added)

        do k=1,num_from_this_row
            edges(1,num_added+k) = i
            edges(2,num_added+k) = g%lists(i)%get_entry(cursor%indx+k)
        enddo

        !! Check that this is right
        !cursor%indx = mod(cursor%indx+num_from_this_row,g%lists(i)%length)

        ! If we returned all nodes from this row, increment the row
        if (num_from_this_row == g%lists(i)%length-cursor%indx) then
            i = i+1
            cursor%indx = 0
        else
            cursor%indx = cursor%indx+num_from_this_row
        endif

        ! Increase the number of edges added
        num_added = num_added+num_from_this_row
    enddo

    cursor%current = cursor%current+num_returned
    cursor%edge(1) = i

end function ll_get_edges



!--------------------------------------------------------------------------!
subroutine ll_add_edge(g,i,j)                                              !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    integer :: cap

    if (.not.g%connected(i,j)) then
        ! Record the old capacity of the list of i's neighbors
        cap = g%lists(i)%capacity

        ! Push the new neighbor node j onto the list of i's neighbors
        call g%lists(i)%push(j)

        ! Change the graph's max degree if need be
        g%max_degree = max(g%max_degree,g%lists(i)%length)

        ! Increment the number of edges and graph capacity
        g%ne = g%ne+1
        g%capacity = g%capacity+g%lists(i)%capacity-cap
    endif

end subroutine ll_add_edge



!--------------------------------------------------------------------------!
subroutine ll_delete_edge(g,i,j)                                           !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    integer :: k,jt,degree,cap

    if (g%connected(i,j)) then
        ! Record the degree of vertex i and capacity
        degree = g%lists(i)%length
        cap = g%lists(i)%capacity

        ! Pop from the list of i's neighbors
        jt = g%lists(i)%pop()

        ! If the vertex jt popped from i's neighbors is not vertex j,
        if (jt/=j) then
            ! find where vertex j was stored and put jt there.
            do k=1,g%lists(i)%length
                if (g%lists(i)%get_entry(k)==j) then
                    call g%lists(i)%set_entry(k,jt)
                endif
            enddo
        endif

        ! If the degree of vertex i was the max degree of the graph, we
        ! need to check that the max degree of the graph hasn't decreased.
        !!Make this a guaranteed O(1) operation somehow
        if (degree==g%max_degree) then
            g%max_degree = 0
            do k=1,g%n
                g%max_degree = max(g%max_degree,g%lists(k)%length)
            enddo
        endif

        ! Decrement the number of edges and capacity of the graph
        g%ne = g%ne-1
        g%capacity = g%capacity+g%lists(i)%capacity-cap
    endif

end subroutine ll_delete_edge



!--------------------------------------------------------------------------!
subroutine ll_graph_left_permute(g,p)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,j,k
    type(dynamic_array) :: lists(g%n)

    do i=1,g%n
        call lists(i)%init(capacity=g%lists(i)%capacity,min_capacity=2)
        do k=1,g%lists(i)%length
            j = g%lists(i)%get_entry(k)
            call lists(i)%push(j)
        enddo
        call g%lists(i)%free()
    enddo

    do i=1,g%n
        call g%lists(p(i))%init(capacity=lists(i)%capacity,min_capacity=2)
        do k=1,lists(i)%length
            j = lists(i)%get_entry(k)
            call g%lists(p(i))%push(j)
        enddo
        call lists(i)%free()
    enddo

end subroutine ll_graph_left_permute



!--------------------------------------------------------------------------!
subroutine ll_graph_right_permute(g,p)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,j,k

    do i=1,g%n
        do k=1,g%lists(i)%length
            j = g%lists(i)%get_entry(k)
            if (j/=0) call g%lists(i)%set_entry(k,p(j))
        enddo
    enddo

end subroutine ll_graph_right_permute



!--------------------------------------------------------------------------!
subroutine ll_left_perm_edge_reorder(g,p,edge_p)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out) :: edge_p(:,:)
    ! local variables
    integer :: i

    ! Initialize the list describing the permutation of the edges
    allocate(edge_p(3,g%n))
    edge_p(1,1) = 1
    do i=1,g%n-1
        edge_p(1,i+1) = edge_p(1,i)+g%lists(i)%length
    enddo

    ! Permute the edges
    call g%left_permute(p)

    ! Finish the list describing the edge permutation
    edge_p(3,1) = 1
    do i=1,g%n-1
        edge_p(3,i+1) = edge_p(3,i)+g%lists(i)%length
    enddo

    do i=1,g%n
        edge_p(2,i) = edge_p(3,p(i))
    enddo

    do i=1,g%n
        edge_p(3,i) = g%lists(p(i))%length
    enddo

end subroutine ll_left_perm_edge_reorder



!--------------------------------------------------------------------------!
subroutine ll_right_perm_edge_reorder(g,p,edge_p)                          !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out) :: edge_p(:,:)

    call g%right_permute(p)
    allocate(edge_p(0,0))

end subroutine ll_right_perm_edge_reorder



!--------------------------------------------------------------------------!
subroutine ll_free(g)                                                      !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer :: i

    do i=1,g%n
        call g%lists(i)%free()
    enddo

    deallocate(g%lists)

    g%n = 0
    g%m = 0
    g%ne = 0
    g%max_degree = 0

end subroutine ll_free



!--------------------------------------------------------------------------!
subroutine ll_dump_edges(g,edges)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: i,k,next

    next = 0
    do i=1,g%n
        do k=1,g%lists(i)%length
            next = next+1
            edges(1,next) = i
            edges(2,next) = g%lists(i)%get_entry(k)
        enddo
    enddo

end subroutine ll_dump_edges






end module ll_graphs

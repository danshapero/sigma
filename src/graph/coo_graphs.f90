module coo_graphs

use graphs
use types, only: dynamic_array

implicit none



!--------------------------------------------------------------------------!
type, extends(graph) :: coo_graph                                          !
!--------------------------------------------------------------------------!
    type(dynamic_array) :: edges(2)
    integer, allocatable, private :: degrees(:)
contains
    !--------------
    ! Constructors
    procedure :: init_const_degree => coo_init_const_degree
    procedure :: init_variable_degree => coo_init_variable_degree
    procedure :: copy => coo_graph_copy

    !-----------
    ! Accessors
    procedure :: degree => coo_degree
    procedure :: neighbors => coo_neighbors
    procedure :: connected => coo_connected
    procedure :: find_edge => coo_find_edge

    !---------------
    ! Edge iterator
    procedure :: make_cursor => coo_make_cursor
    procedure :: get_edges => coo_get_edges

    !----------
    ! Mutators
    procedure :: add_edge => coo_add_edge
    procedure :: delete_edge => coo_delete_edge
    procedure :: left_permute => coo_graph_left_permute
    procedure :: right_permute => coo_graph_right_permute
    procedure :: compress => coo_graph_compress
    procedure :: decompress => coo_graph_decompress

    !-------------
    ! Destructors
    procedure :: free => coo_free

    !--------------------------
    ! Testing, debugging & I/O
    procedure :: dump_edges => coo_dump_edges
end type coo_graph





contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine coo_init_const_degree(g,n,m,degree)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, degree
    ! local variables
    integer :: ne

    g%n = n
    allocate(g%degrees(n))
    g%degrees = 0

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    if (present(degree)) then
        ne = degree*g%n
    else
        ne = max(g%m,g%n)
    endif

    call g%edges(1)%init(capacity=ne)
    call g%edges(2)%init(capacity=ne)

    ! At initialization, the number of edges and max degree is zero
    g%ne = 0
    g%max_degree = 0
    g%capacity = ne

    ! Mark the graph as mutable
    g%mutable = .true.

end subroutine coo_init_const_degree



!--------------------------------------------------------------------------!
subroutine coo_init_variable_degree(g,n,m,degrees)                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: n, degrees(:)
    integer, intent(in), optional :: m
    ! local variables
    integer :: ne

    g%n = n
    allocate(g%degrees(n))
    g%degrees = 0

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    ! The sum of the degree list gives a bound on how much space to allocate
    ne = sum(degrees)

    call g%edges(1)%init(capacity=ne)
    call g%edges(2)%init(capacity=ne)

    ! At initialization, the number of edges and max degree is zero
    g%ne = 0
    g%max_degree = 0
    g%capacity = ne

    ! Mark the graph as mutable
    g%mutable = .true.

end subroutine coo_init_variable_degree



!--------------------------------------------------------------------------!
subroutine coo_graph_copy(g,h)                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    class(graph), intent(in)        :: h
    ! local variables
    integer :: i,j,k,n,num_blocks,num_returned,edges(2,64)
    type(graph_edge_cursor) :: cursor

    ! Mark the graph as mutable
    g%mutable = .true.

    ! Copy all the attributes of h to g
    g%n = h%n
    g%m = h%m
    g%ne = h%ne
    g%max_degree = h%max_degree

    g%capacity = g%ne

    allocate(g%degrees(g%capacity))

    ! Allocate space in the two dynamic arrays for the edges of g
    call g%edges(1)%init(capacity=g%capacity)
    call g%edges(2)%init(capacity=g%capacity)

    ! Make an edge iterator for the copied graph h
    cursor = h%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    ! Iterate through all the edges of h
    do n=1,num_blocks
        ! Get a chunk of edges of h
        edges = h%get_edges(cursor,64,num_returned)

        ! Add each edge from the chunk into g
        do k=1,num_returned
            i = edges(1,k)
            j = edges(2,k)

            if (i/=0 .and. j/=0) then
                call g%edges(1)%push(i)
                call g%edges(2)%push(j)

                g%degrees(i) = g%degrees(i)+1
            endif
        enddo
    enddo

end subroutine coo_graph_copy




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function coo_degree(g,i) result(d)                                         !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: d

    d = g%degrees(i)

end function coo_degree



!--------------------------------------------------------------------------!
subroutine coo_neighbors(g,i,nbrs)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer, intent(in) :: i
    integer, intent(out) :: nbrs(:)
    ! local variables
    integer :: k,next

    nbrs = 0
    next = 0
    do k=1,g%ne
        if (g%edges(1)%get_entry(k)==i) then
            next = next+1
            nbrs(next) = g%edges(2)%get_entry(k)
        endif
    enddo

end subroutine coo_neighbors



!--------------------------------------------------------------------------!
function coo_connected(g,i,j)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: coo_connected
    ! local variables
    integer :: k

    coo_connected = .false.

    do k=1,g%ne
        if (g%edges(1)%get_entry(k)==i.and.g%edges(2)%get_entry(k)==j) then
            coo_connected = .true.
        endif
    enddo

end function coo_connected



!--------------------------------------------------------------------------!
function coo_find_edge(g,i,j)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer, intent(in) :: i,j
    integer :: coo_find_edge
    ! local variables
    integer :: k

    coo_find_edge = -1

    do k=1,g%ne
        if (g%edges(1)%get_entry(k)==i.and.g%edges(2)%get_entry(k)==j) then
            coo_find_edge = k
        endif
    enddo

end function coo_find_edge




!==========================================================================!
!==== Edge iterator                                                    ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function coo_make_cursor(g,thread) result(cursor)                          !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(in) :: g
    integer, intent(in) :: thread
    type(graph_edge_cursor) :: cursor

    cursor%start = 1
    cursor%final = g%ne
    cursor%current = 0
    cursor%edge = [g%edges(1)%get_entry(1), g%edges(2)%get_entry(1)]

end function coo_make_cursor



!--------------------------------------------------------------------------!
function coo_get_edges(g,cursor,num_edges,num_returned) result(edges)      !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(in) :: g
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    integer :: edges(2,num_edges)

    ! Set up the returned edges to be 0
    edges = 0

    ! Count how many edges we're actually going to return; we'll either
    ! return how many edges the user asked for, or, if that amount would
    ! go beyond the final edge that the cursor is allowed to access, the
    ! all of the remaining edges
    num_returned = min(num_edges,cursor%final-cursor%current)

    ! Fill the edges array with the right slice from the graph's edges
    edges(1,1:num_returned) = &
        & g%edges(1)%array(cursor%current+1:cursor%current+num_returned)
    edges(2,1:num_returned) = &
        & g%edges(2)%array(cursor%current+1:cursor%current+num_returned)

    ! Move the cursor's current edge ahead to the last one we returned
    cursor%current = cursor%current+num_returned

end function coo_get_edges




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine coo_add_edge(g,i,j)                                             !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: i,j

    if (.not.g%mutable) then
        print *, 'Attempted to add an edge to an immutable COO graph'
        print *, 'Terminating.'
        call exit(1)
    endif

    if (.not.g%connected(i,j)) then
        ! Add in the new edge
        call g%edges(1)%push(i)
        call g%edges(2)%push(j)

        ! Increase the number of edges and the graph capacity
        g%ne = g%ne+1
        g%capacity = g%edges(1)%capacity

        ! If the degree of node i is now the greatest of all nodes in
        ! the graph, update the degree accordingly
        g%degrees(i) = g%degrees(i)+1
        g%max_degree = max(g%max_degree,g%degrees(i))
    endif

end subroutine coo_add_edge



!--------------------------------------------------------------------------!
subroutine coo_delete_edge(g,i,j)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k,it,jt

    if (.not.g%mutable) then
        print *, 'Attempted to delete an edge from an immutable COO graph'
        print *, 'Terminating.'
        call exit(1)
    endif

    k = g%find_edge(i,j)

    if (k/=-1) then
        it = g%edges(1)%pop()
        jt = g%edges(2)%pop()

        if (k<g%ne .and. g%ne>1) then
            call g%edges(1)%set_entry(k,it)
            call g%edges(2)%set_entry(k,jt)
        endif

        ! Decrement the number of edges and the graph capacity
        g%ne = g%ne-1
        g%capacity = g%edges(1)%capacity

        ! Change the degree of vertex i
        g%degrees(i) = g%degrees(i)-1

        ! If vertex i had the greatest degree in the graph, check to
        ! make sure the max degree hasn't decreased
        if (g%degrees(i)+1==g%max_degree) then
            g%max_degree = maxval(g%degrees)
        endif
    endif

end subroutine coo_delete_edge



!--------------------------------------------------------------------------!
subroutine coo_graph_left_permute(g,p,edge_p)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i,k

    do k=1,g%ne
        i = g%edges(1)%get_entry(k)
        if (i/=0) call g%edges(1)%set_entry(k,p(i))
    enddo

    g%degrees = 0
    do k=1,g%ne
        i = g%edges(1)%get_entry(k)
        g%degrees(i) = g%degrees(i)+1
    enddo

    if (present(edge_p)) allocate(edge_p(0,0))

end subroutine coo_graph_left_permute



!--------------------------------------------------------------------------!
subroutine coo_graph_right_permute(g,p,edge_p)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i,k

    do k=1,g%ne
        i = g%edges(2)%get_entry(k)
        if (i/=0) call g%edges(2)%set_entry(k,p(i))
    enddo

    if (present(edge_p)) allocate(edge_p(0,0))

end subroutine coo_graph_right_permute



!--------------------------------------------------------------------------!
subroutine coo_graph_compress(g,edge_p)                                    !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(inout) :: g
    integer, allocatable, intent(inout), optional :: edge_p(:,:)

    if (present(edge_p)) allocate(edge_p(0,0))

    ! COO graphs cannot have their storage compressed.

    ! Mark the graph as immutable.
    g%mutable = .false.

end subroutine coo_graph_compress



!--------------------------------------------------------------------------!
subroutine coo_graph_decompress(g)                                         !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(inout) :: g

    g%mutable = .true.

end subroutine coo_graph_decompress



!--------------------------------------------------------------------------!
subroutine coo_free(g)                                                     !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(inout) :: g

    deallocate(g%edges(1)%array,g%edges(2)%array,g%degrees)

    g%n = 0
    g%m = 0
    g%ne = 0
    g%capacity = 0
    g%max_degree = 0

    g%mutable = .true.

end subroutine coo_free



!--------------------------------------------------------------------------!
subroutine coo_dump_edges(g,edges)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: k

    do k=1,g%ne
        edges(1,k) = g%edges(1)%get_entry(k)
        edges(2,k) = g%edges(2)%get_entry(k)
    enddo

end subroutine coo_dump_edges



end module coo_graphs

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
    procedure :: add_edge => ellpack_add_edge
    procedure :: delete_edge => ellpack_delete_edge
    procedure :: left_permute => ellpack_graph_left_permute, &
                & right_permute => ellpack_graph_right_permute
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

    allocate( g%node(g%max_degree,ne) )

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

    nbrs = g%node(:,i)

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

    do k=1,g%max_degree
        if (g%node(k,i)==j) ellpack_find_edge = k
    enddo

end function ellpack_find_edge



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
    integer :: i,k,node(g%max_degree,g%n)

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

    g%node(i,g%max_degree+1) = j

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

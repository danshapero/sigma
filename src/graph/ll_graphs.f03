module ll_graphs

use abstract_graphs
use types

implicit none


!--------------------------------------------------------------------------!
type, extends(graph) :: ll_graph                                           !
!--------------------------------------------------------------------------!
    type(linked_list), allocatable :: lists(:)
contains
    procedure :: init => ll_init
    procedure :: neighbors => ll_neighbors
    procedure :: connected => ll_connected
    procedure :: find_edge => ll_find_edge
    procedure :: add_edge => ll_add_edge
    procedure :: delete_edge => ll_delete_edge
    procedure :: free => ll_free
    procedure :: dump_edges => ll_dump_edges

end type ll_graph





contains



!--------------------------------------------------------------------------!
subroutine ll_init(g,n,m,edges)                                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, edges(:,:)
    ! local variables
    integer :: k,ne

    g%n = n
    allocate(g%lists(n))

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    if (present(edges)) then
        ne = size(edges,2)
        g%ne = ne

        do k=1,ne
            call g%lists(edges(1,k))%append(edges(2,k))
        enddo

        g%max_degree = 0
        do k=1,n
            g%max_degree = max(g%max_degree,g%lists(k)%length)
        enddo
    else
        ne = 0
        g%ne = ne
        g%max_degree = 0
    endif

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
        nbrs(k) = g%lists(i)%get_value(k)
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
        if (g%lists(i)%get_value(k)==j) then
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
    integer :: k

    ll_find_edge = -1
    do k=1,g%lists(i)%length
        if (g%lists(i)%get_value(k)==j) then
            ll_find_edge = k
            exit
        endif
    enddo

end function ll_find_edge



!--------------------------------------------------------------------------!
subroutine ll_add_edge(g,i,j)                                              !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: i,j

    if (.not.g%connected(i,j)) then
        call g%lists(i)%prepend(j)
        if (g%lists(i)%length>g%max_degree) then
            g%max_degree = g%lists(i)%length
        endif
        g%ne = g%ne+1
    endif

end subroutine ll_add_edge



!--------------------------------------------------------------------------!
subroutine ll_delete_edge(g,i,j)                                           !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    integer :: k,degree

    degree = g%lists(i)%length
    call g%lists(i)%delete_value(j)
    if (g%lists(i)%length<degree) then
        g%ne = g%ne-1
    endif

    if (degree==g%max_degree) then
        g%max_degree = 0
        do k=1,g%n
            g%max_degree = max(g%max_degree,g%lists(k)%length)
        enddo
    endif

end subroutine ll_delete_edge



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
            edges(2,next) = g%lists(i)%get_value(k)
        enddo
    enddo

end subroutine ll_dump_edges






end module ll_graphs

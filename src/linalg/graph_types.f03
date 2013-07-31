module graph_types

use graphs

implicit none



!--------------------------------------------------------------------------!
type, extends(graph) :: csr_graph                                          !
!--------------------------------------------------------------------------!
    integer, allocatable :: ia(:), ja(:)
contains
    procedure :: init => csr_init
    procedure :: neighbors => csr_neighbors
    procedure :: connected => csr_connected
    procedure :: find_edge => csr_find_edge
    procedure :: add_edge  => csr_add_edge
    procedure :: delete_edge => csr_delete_edge
    ! auxiliary routines
    !procedure :: sort_ja
end type csr_graph





contains



!--------------------------------------------------------------------------!
subroutine csr_init(g,n,m,edges)                                           !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, edges(:,:)
    ! local variables
    integer :: i,k,work(n),ne

    g%n = n
    allocate( g%ia(n+1) )

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

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
    endif

end subroutine csr_init



!--------------------------------------------------------------------------!
subroutine csr_neighbors(g,i,nbrs)                                         !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_graph), intent(in) :: g
    integer, intent(in) :: i
    integer, intent(out) :: nbrs(:)
    ! local variables
    integer :: start,finish,degree

    degree = g%ia(i+1)-g%ia(i)
    start = g%ia(i)
    finish = g%ia(i+1)-1

    nbrs = 0
    nbrs(1:degree) = g%ja(start:finish)

end subroutine csr_neighbors



!--------------------------------------------------------------------------!
function csr_connected(g,i,j)                                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: csr_connected
    ! local variables
    integer :: k

    csr_connected = .false.

    do k=g%ia(i),g%ia(i+1)-1
        if (g%ja(k)==j) csr_connected = .true.
    enddo

end function csr_connected



!--------------------------------------------------------------------------!
function csr_find_edge(g,i,j)                                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_graph), intent(in) :: g
    integer, intent(in) :: i,j
    integer :: csr_find_edge
    ! local variables
    integer :: k

    csr_find_edge = g%ia(i)

    do k=g%ia(i),g%ia(i+1)-1
        if (g%ja(k)<j) then
            csr_find_edge = k+1
        endif
    enddo

end function csr_find_edge



!--------------------------------------------------------------------------!
subroutine csr_add_edge(g,i,j)                                             !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k,indx,ja_temp(g%ne)

    if (.not.g%connected(i,j)) then
        ! Find the index in the list of edges where the new edge is to be
        ! inserted
        indx = g%find_edge(i,j)

        ja_temp = g%ja

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
    endif

end subroutine csr_add_edge



!--------------------------------------------------------------------------!
subroutine csr_delete_edge(g,i,j)                                          !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(csr_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k,indx,ja_temp(g%ne)

    if (g%connected(i,j)) then
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
    endif

end subroutine csr_delete_edge





end module graph_types

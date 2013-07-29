module graphs

implicit none



!--------------------------------------------------------------------------!
type :: graph                                                              !
!--------------------------------------------------------------------------!
    ! Data for graph storage
    integer :: nrow,ncol,ne,max_degree
    integer, allocatable :: ia(:),ja(:)
    logical :: directed
    ! Object-bound procedures
    procedure(init_ifc), pointer        :: init
    procedure(neighbors_ifc), pointer   :: neighbors
!    procedure(add_edge_ifc), pointer    :: add_edge
!    procedure(delete_edge_ifc), pointer :: delete_edge
    procedure(connected_ifc), pointer   :: connected
end type graph


!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine init_ifc(g,edges,nrow,ncol)
        import :: graph
        class(graph), intent(inout) :: g
        integer, intent(in) :: edges(:,:), nrow, ncol
    end subroutine init_ifc

    function neighbors_ifc(g,i)
        import :: graph
        class(graph), intent(in) :: g
        integer, intent(in) :: i
        integer :: neighbors_ifc(g%max_degree)
    end function neighbors_ifc

    function connected_ifc(g,i,j)
        import :: graph
        class(graph), intent(in) :: g
        integer, intent(in) :: i,j
        logical :: connected_ifc
    end function connected_ifc


end interface



contains



function get_graph(storage) result(g)
    implicit none
    character(len=*), intent(in), optional :: storage
    type(graph) :: g

    g%init          => coo_init
    g%neighbors     => coo_neighbors
    g%connected     => coo_connected

end function get_graph




!--------------------------------------------------------------------------!
! Graph procedures for the coordinate format                               !
!--------------------------------------------------------------------------!

subroutine coo_init(g,edges,nrow,ncol)
    implicit none
    ! input/output variables
    class(graph), intent(inout) :: g
    integer, intent(in) :: edges(:,:), nrow, ncol
    ! local variables
    integer :: i,n,degree(nrow)

    g%nrow = nrow
    g%ncol = ncol
    g%ne = size(edges,2)

    allocate(g%ia(g%ne),g%ja(g%ne))

    g%ia = edges(1,:)
    g%ja = edges(2,:)

    degree = 0
    do n=1,g%ne
        i = g%ia(n)
        degree(i) = degree(i)+1
    enddo
    g%max_degree = maxval(degree)
    
end subroutine coo_init



function coo_neighbors(g,i)
    implicit none
    ! input/output variables
    class(graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: coo_neighbors(g%max_degree)
    ! local variables
    integer :: k,next

    coo_neighbors = 0
    next = 0
    do k=1,g%ne
        if (g%ia(k)==i) then
            next = next+1
            coo_neighbors(next) = g%ja(k)
        endif
    enddo

end function coo_neighbors



function coo_connected(g,i,j)
    implicit none
    ! input/output variables
    class(graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: coo_connected
    ! local variables
    integer :: k

    coo_connected = .false.

    do k=1,g%ne
        if ( g%ia(k)==i .and. g%ja(k)==j ) coo_connected = .true.
    enddo

end function coo_connected



end module graphs

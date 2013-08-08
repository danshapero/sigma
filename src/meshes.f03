module meshes

use graphs
use ll_graphs
use cs_graphs
use conversions

implicit none



type :: mesh
    integer :: d, n_nodes, n_edges, n_faces, n_vols, &
        & nodes, edges, faces, elements
    real(dp), allocatable :: x(:,:)
    integer, allocatable :: bnd(:), ebnd(:), fbnd(:)
    type(graph_pointer), allocatable :: gp(:,:)
end type



contains


!==========================================================================!
!==========================================================================!
!==== Routines for reading meshes in the Triangle/Tetgen format        ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine read_triangle_mesh(g,x,bnd,ele,filename)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout), allocatable :: g
    real(dp), intent(inout), allocatable :: x(:,:)
    integer, intent(inout), allocatable :: bnd(:), ele(:,:)
    character(len=*), intent(in) :: filename
    ! local variables
    integer :: i,j,k,n,n_node,n_ele,d
    logical :: file_exists

    inquire(file=trim(filename)//'.node', exist=file_exists)
    if (file_exists) then
        call read_triangle_nodes(x,bnd,filename)
    else
        print *, trim(filename)//'.node not found'
        call exit(1)
    endif
    d = size(x,1)
    n_node = size(x,2)

    inquire(file=trim(filename)//'.ele', exist=file_exists)
    if (file_exists) then
        call read_triangle_elements(ele,filename)
    else
        print *, trim(filename)//'.ele not found'
        call exit(1)
    endif
    n_ele = size(ele,2)

    allocate(ll_graph::g)
    call g%init(n_node)

    do n=1,n_ele
        do k=1,d+1
            i = ele(n,k)
            j = ele(n,mod(k,d+1)+1)
            call g%add_edge(i,j)
            call g%add_edge(j,i)
        enddo
    enddo
    do i=1,n_node
        call g%add_edge(i,i)
    enddo

    call convert(g)

end subroutine read_triangle_mesh



!--------------------------------------------------------------------------!
subroutine read_triangle_nodes(x,bnd,filename)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    real(dp), intent(inout), allocatable :: x(:,:)
    integer, intent(inout), allocatable :: bnd(:)
    character(len=*), intent(in) :: filename
    ! local variables
    logical :: file_exists
    integer :: dummy, n, n_nodes, d

    inquire(file=trim(filename)//'.node', exist=file_exists)
    if (file_exists) then
        open(unit=10, file=trim(filename)//'.node')
        read(10,*) n_nodes, d
        allocate(x(d,n_nodes), bnd(n_nodes))

        do n=1,n_nodes
            read(10,*) dummy, x(:,n), bnd(n)
        enddo

        close(10)
    else
        print *, trim(filename)//'.node not found! Exiting'
        call exit(1)
    endif

end subroutine read_triangle_nodes



!--------------------------------------------------------------------------!
subroutine read_triangle_edges(edges,ebnd,filename)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    integer, intent(inout), allocatable :: edges(:,:), ebnd(:)
    character(len=*), intent(in) :: filename
    ! local variables
    logical :: file_exists
    integer :: dummy, n, n_edges

    inquire(file=trim(filename)//'.edge', exist=file_exists)
    if (file_exists) then
        open(unit=10, file=trim(filename)//'.edge')

        read(10,*) n_edges
        allocate(edges(2,n_edges),ebnd(n_edges))

        do n=1,n_edges
            read(10,*) dummy, edges(1:2,n), ebnd(n)
        enddo

        close(10)
    else
        print *, trim(filename)//'.edge not found! Exiting'
        call exit(1)
    endif

end subroutine read_triangle_edges



!--------------------------------------------------------------------------!
subroutine read_triangle_faces(faces,fbnd,filename)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    integer, intent(inout), allocatable :: faces(:,:), fbnd(:)
    character(len=*), intent(in) :: filename
    ! local variables
    logical :: file_exists
    integer :: dummy, n, n_faces

    inquire(file=trim(filename)//'.face', exist=file_exists)
    if (file_exists) then
        open(unit=10, file=trim(filename)//'.face')

        read(10,*) n_faces
        allocate(faces(3,n_faces),fbnd(n_faces))

        do n=1,n_faces
            read(10,*) dummy, faces(1:3,n), fbnd(n)
        enddo

        close(10)
    else
        print *, trim(filename)//'.face not found! Exiting'
        call exit(1)
    endif

end subroutine read_triangle_faces



!--------------------------------------------------------------------------!
subroutine read_triangle_elements(elements,filename)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    integer, intent(inout), allocatable :: elements(:,:)
    character(len=*), intent(in) :: filename
    ! local variables
    logical :: file_exists
    integer :: dummy, d, n, n_ele

    inquire(file=trim(filename)//'.ele', exist=file_exists)
    if (file_exists) then
        open(unit=10, file=trim(filename)//'.ele')

        read(10,*) n_ele, d
        allocate(elements(d,n_ele))

        do n=1,n_ele
            read(10,*) dummy, elements(1:d,n)
        enddo

        close(10)
    else
        print *, trim(filename)//'.ele not found! Exiting'
        call exit(1)
    endif

end subroutine read_triangle_elements



!--------------------------------------------------------------------------!
subroutine read_triangle_neighbors(neighbors,filename)                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    integer, intent(inout), allocatable :: neighbors(:,:)
    character(len=*), intent(in) :: filename
    ! local variables
    logical :: file_exists
    integer :: dummy, d, n, n_ele

    inquire(file=trim(filename)//'.neigh', exist=file_exists)
    if (file_exists) then
        open(unit=10, file=trim(filename)//'.neigh')

        read(10,*) n_ele, d
        allocate(neighbors(d,n_ele))

        do n=1,n_ele
            read(10,*) dummy, neighbors(1:d,n)
        enddo

        close(10)
    else
        print *, trim(filename)//'.neigh not found! Exiting'
        call exit(1)
    endif


end subroutine read_triangle_neighbors







!==========================================================================!
!==========================================================================!
!==== Routines for reading meshes in gmsh format                       ====!
!==========================================================================!
!==========================================================================!







!==========================================================================!
!==========================================================================!
!==== Routines for reading meshes in CGAL format                       ====!
!==========================================================================!
!==========================================================================!





end module meshes

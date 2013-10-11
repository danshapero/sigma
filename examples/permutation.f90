program permutation

use fempack
use fem

implicit none

    character(len=64) :: mesh_filename, out_filename

    ! variables for constructing the computational mesh
    class(graph), pointer :: g
    real(dp), allocatable :: x(:,:)
    integer, allocatable :: ele(:,:), bnd(:)

    integer :: nn, ne, n
    integer, allocatable :: p(:), edges(:,:)


!--------------------------------------------------------------------------!
! Read in the mesh from a file and assemble the graph                      !
!--------------------------------------------------------------------------!
    call getarg(1,mesh_filename)
    call getarg(2,out_filename)

    call read_triangle_mesh(g,x,bnd,ele,trim(mesh_filename))
    nn = g%n
    ne = g%ne

    allocate(p(nn))


!--------------------------------------------------------------------------!
! Compute the breadth-first reordering of the graph                        !
!--------------------------------------------------------------------------!
    call breadth_first_search(g,p)

    call g%right_permute(p)
    call g%left_permute(p)


!--------------------------------------------------------------------------!
! Write the result to a file                                               !
!--------------------------------------------------------------------------!
    allocate(edges(2,ne))
    call g%dump_edges(edges)

    open(unit=10,file=trim(out_filename))
    write(10,*) nn, ne
    do n=1,ne
        write(10,*) edges(:,n)
    enddo
    close(10)


end program permutation

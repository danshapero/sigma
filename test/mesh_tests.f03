program mesh_tests

use fempack

implicit none


    class(graph), pointer :: g
    real(dp), allocatable :: x(:,:)
    integer, allocatable :: bnd(:), ele(:,:)

    call read_triangle_mesh(g,x,bnd,ele,'examples/meshes/circle.1')


end program mesh_tests

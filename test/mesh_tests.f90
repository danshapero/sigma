program mesh_tests

use fempack

implicit none


    class(graph), pointer :: g
    real(dp), allocatable :: x(:,:)
    integer, allocatable :: bnd(:), ele(:,:)
    character(len=32) :: directory

    call get_environment_variable("FEMPACK",directory)

    if (trim(directory)=="") then
        print *, "Need to have the environment variable FEMPACK "
        print *, "set to the fempack source directory."
        call exit(1)
    endif

    call read_triangle_mesh(g,x,bnd,ele,trim(directory)// &
        & '/examples/meshes/circle.1')


end program mesh_tests

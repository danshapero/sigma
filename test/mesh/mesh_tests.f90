program mesh_tests

    use meshes

    implicit none

    ! computational mesh
    type(tri_mesh) :: mesh

    ! some other locals
    integer :: i,j


!--------------------------------------------------------------------------!
! Test reading the mesh                                                    !
!--------------------------------------------------------------------------!
    call read_mesh("../../meshes/circle.1",mesh)


end program mesh_tests

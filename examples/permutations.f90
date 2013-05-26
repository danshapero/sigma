program permutations

    use meshes
    use linalg
    use fem
    use netcdf

    implicit none

    ! computational mesh
    type (tri_mesh) :: mesh

    ! stiffness matrix
    class (sparse_matrix), allocatable :: A

    ! permutation
    integer, allocatable :: p(:)

    ! some other locals
    integer :: maxcolor


!--------------------------------------------------------------------------!
! Load the mesh and construct the system matrix                            !
!--------------------------------------------------------------------------!
    call read_mesh("../meshes/circle.1              ",mesh)

    allocate(csr_matrix::A)

    call assemble(mesh,A)
    call stiffness_matrix(mesh,A,1.d0)

    call A%write_to_file("a")


!--------------------------------------------------------------------------!
! Breadth-first search permutation                                         !
!--------------------------------------------------------------------------!
    allocate( p(A%nrow) )
    call bfs(A,p)
    call A%permute(p)
    call A%write_to_file("a_bfs")



!--------------------------------------------------------------------------!
! Multicolor permutation                                                   !
!--------------------------------------------------------------------------!
    p = 0
    call greedy_multicolor(A,p,maxcolor)
    call A%permute(p)
    call A%write_to_file("a_multicolor")



end program permutations

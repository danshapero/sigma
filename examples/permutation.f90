program permutation

use fempack
use fem

implicit none

    character(len=64) :: mesh_filename, out_filename

    ! variables for constructing the computational mesh
    class(graph), pointer :: g
    real(dp), allocatable :: x(:,:)
    integer, allocatable :: ele(:,:), bnd(:)

    integer :: nn, ne, k, n, num_colors, min_nbr, max_nbr, bandwidth
    integer, allocatable :: p(:), colors(:), ptrs(:), nbrs(:)


!--------------------------------------------------------------------------!
! Read in the mesh from a file and assemble the graph                      !
!--------------------------------------------------------------------------!
    call getarg(1,mesh_filename)
    call getarg(2,out_filename)

    call read_triangle_mesh(g,x,bnd,ele,trim(mesh_filename))
    nn = g%n
    ne = g%ne

    allocate(p(nn), nbrs(g%max_degree))


!--------------------------------------------------------------------------!
! Compute the breadth-first reordering of the graph                        !
!--------------------------------------------------------------------------!
    call breadth_first_search(g,p)

    call g%right_permute(p)
    call g%left_permute(p)

    bandwidth = 0
    do n=1,g%n
        call g%neighbors(n,nbrs)

        min_nbr = g%n+1
        max_nbr = 0
        do k=1,g%max_degree
            max_nbr = max(max_nbr,nbrs(k))
            min_nbr = min(min_nbr,nbrs(k))
        enddo

        bandwidth = max(bandwidth,(max_nbr-min_nbr)/2)
    enddo

    print *, float(bandwidth)/g%max_degree

    call g%write_to_file(trim(out_filename)//'_bfs.txt')



!--------------------------------------------------------------------------!
! Find the multicolor ordering of the graph                                !
!--------------------------------------------------------------------------!
    allocate(ptrs(g%max_degree+1))
    call greedy_color_ordering(g,p,ptrs,num_colors)

    call g%right_permute(p)
    call g%left_permute(p)

    allocate(colors(num_colors))
    do n=1,num_colors
        colors(n) = ptrs(n+1)-ptrs(n)
    enddo
    print *, float(minval(colors))/maxval(colors)

    call g%write_to_file(trim(out_filename)//'_mc.txt')



end program permutation

!==========================================================================!
program matrix_example_2                                                   !
!==========================================================================!
!==== Example program demonstrating eigenvalue approximation of sparse ====!
!==== matrices using the Lanczos process.                              ====!
!==========================================================================!

use sigma

implicit none

    ! a linked-list graph and variables for generating it randomly
    class(graph), pointer :: g
    real(dp) :: z(2), p

    ! a sparse_matrix
    type(sparse_matrix) :: A

    ! some integer indices
    integer :: i,j,k,x,y,nx,ny,d,di,dj,n
    integer, allocatable :: neighbors(:)

    ! some vectors
    real(dp), allocatable :: lambda(:), V(:,:)


    allocate(cs_graph::g)

    n = 128
    nx = 32
    ny = 32

    allocate(lambda(n),V(nx*ny,n))
    lambda = 0
    V = 0


    ! initialize a random seed and choose the probability of two vertices
    ! being connected
    call init_seed()
    p = 0.25



    !----------------------------------------------------------------------!
    ! Set up a graph representing the integer lattice on a grid of         !
    ! nx x ny points.                                                      !
    !----------------------------------------------------------------------!
    call g%init(nx*ny,degree=5)

    ! Make every (x,y) connected to (x,y+1), (x,y-1), (x+1,y) & (x-1,y)
    do x = 1,nx
        do y = 1,ny
            i = nx*(x-1)+y
            call g%add_edge(i,i)

            ! Add edge from (x,y) to (x,y+1) and vice versa
            j = nx*(x-1)+mod(y,ny)+1
            call g%add_edge(i,j)
            call g%add_edge(j,i)

            ! (x,y) to (x+1,y) and v.v.
            j = nx*mod(x,nx)+y
            call g%add_edge(i,j)
            call g%add_edge(j,i)
        enddo
    enddo

    write(*,10) g%ne/2
10  format('Done building initial graph;  ',i8,' edges.')


    !----------------------------------------------------------------------!
    ! Randomly delete edges.                                               !
    !----------------------------------------------------------------------!
    do x=1,nx
        do y=1,ny
            ! Harvest two random numbers
            call random_number(z)

            ! Find the index for the point (x,y)
            i = nx*(x-1)+y

            ! If the first number is less than p, delete the connection
            ! to (x,y+1)
            j = nx*(x-1)+mod(y,ny)+1
            if (z(1)<p) then
                call g%delete_edge(i,j)
                call g%delete_edge(j,i)
            endif

            ! If the second number is less than p, delete the connection
            ! to (x+1,y)
            j = nx*mod(x,nx)+y
            if (z(2)<p) then
                call g%delete_edge(i,j)
                call g%delete_edge(j,i)
            endif
        enddo
    enddo

    write(*,20) g%ne/2
20  format('Done randomly removing edges; ',i8,' edges remain.')



    !----------------------------------------------------------------------!
    ! Compress the graph storage.                                          !
    !----------------------------------------------------------------------!
    call g%compress()



    !----------------------------------------------------------------------!
    ! Make the Laplacian matrix of g                                       !
    !----------------------------------------------------------------------!
    call A%init(nx*ny,nx*ny,'row',g)
    allocate(neighbors(g%max_degree))
    do i=1,A%nrow
        di = g%degree(i)-1
        call A%set_value(i,i,1.0_dp)
        call g%neighbors(i,neighbors)
        do k=1,di+1
            j = neighbors(k)
            if (j/=0 .and. j/=i) then
                dj = g%degree(j)-1
                call A%set_value(i,j,-1.0_dp/dsqrt(1.0_dp*di*dj))
            endif
        enddo
    enddo



    !----------------------------------------------------------------------!
    ! Call the eigensolver                                                 !
    !----------------------------------------------------------------------!
    call eigensolve(A,lambda,V,n)
    write(*,30) lambda(2), lambda(3), lambda(4)
30  format('Three smallest eigenvalues of graph: ',F9.6,', ',F9.6,', ',F9.6)
    write(*,40) (lambda(3)-lambda(2))/lambda(2)
40  format('Eigenvalue gap: ',F12.6)

    open(unit=10,file="eigs.txt")
    write(10,*) nx*ny, n
    do i=1,nx*ny
        write(10,*) V(i,:)
    enddo
    close(10)

    open(unit=20,file="graph.txt")
    write(20,*) nx*ny
    do i=1,nx*ny
        call g%neighbors(i,neighbors)
        d = g%degree(i)
        write(20,*) i,neighbors(1:d)
    enddo
    close(20)


end program matrix_example_2

program block_matrix_tests

use fempack

implicit none

    type(block_matrix) :: A
    class(graph), pointer :: g
    class(sparse_matrix), pointer :: M
    integer :: i,j,k,n
    real(dp) :: z
    type(vector) :: x,y
    integer :: n_nodes, n_triangles
    integer, allocatable :: triangles(:,:), neigh(:,:), nbrs(:)

    ! Hard-code the mesh geometry
    n_nodes = 15
    n_triangles = 18
    allocate(triangles(3,18),neigh(3,18))

    triangles(:,1)  = [ 1, 3, 4]
    triangles(:,2)  = [ 1, 2, 4]
    triangles(:,3)  = [ 3, 4, 5]
    triangles(:,4)  = [ 4, 5, 8]
    triangles(:,5)  = [ 4, 7, 8]
    triangles(:,6)  = [ 3, 4, 7]
    triangles(:,7)  = [ 5, 8, 9]
    triangles(:,8)  = [ 5, 6, 9]
    triangles(:,9)  = [ 2, 5, 6]
    triangles(:,10) = [ 6, 9,10]
    triangles(:,11) = [ 9,10,13]
    triangles(:,12) = [ 9,13,14]
    triangles(:,13) = [ 9,14,12]
    triangles(:,14) = [ 8,12, 9]
    triangles(:,15) = [ 7,12, 8]
    triangles(:,16) = [ 7,11,12]
    triangles(:,17) = [11,15,12]
    triangles(:,18) = [15,14,12]

    neigh(:,1)  = [ 2,-1, 6]
    neigh(:,2)  = [ 1, 3,-1]
    neigh(:,3)  = [ 2, 4, 9]
    neigh(:,4)  = [ 3, 5, 7]
    neigh(:,5)  = [ 4, 6,15]
    neigh(:,6)  = [ 5, 1,-1]
    neigh(:,7)  = [ 4,14, 8]
    neigh(:,8)  = [ 7,10, 9]
    neigh(:,9)  = [ 3, 8,-1]
    neigh(:,10) = [ 8,11,-1]
    neigh(:,11) = [10,12,-1]
    neigh(:,12) = [11,13,-1]
    neigh(:,13) = [12,14,18]
    neigh(:,14) = [ 7,15,13]
    neigh(:,15) = [ 5,16,14]
    neigh(:,16) = [17,15,-1]
    neigh(:,17) = [18,16,-1]
    neigh(:,18) = [13,17,-1]


    !----------------------------------------------------------------------!
    ! Set up the block matrix                                              !
    !----------------------------------------------------------------------!
    call A%init(n_nodes+n_triangles,n_nodes+n_triangles,'row')
    call A%assemble([n_nodes,n_triangles],[n_nodes,n_triangles])

    ! Create the graph and matrix for the nodes
    allocate(ll_graph::g)
    call g%init(n_nodes,n_nodes)

    do i=1,n_nodes
        call g%add_edge(i,i)
    enddo

    do n=1,n_triangles
        do k=1,3
            i = triangles(k,n)
            j = triangles(mod(k,3)+1,n)
            call g%add_edge(i,j)
            call g%add_edge(j,i)
        enddo
    enddo

    call convert(g,'cs')

    allocate(cs_matrix::M)
    call M%init(n_nodes,n_nodes,'row',g)

    allocate(nbrs(g%max_degree))
    do i=1,n_nodes
        call g%neighbors(i,nbrs)
        do k=1,g%max_degree
            j = nbrs(k)
            if (j/=0 .and. j/=i) then
                call M%set_value(i,j,-1.0_dp)
                call M%add_value(i,i,1.0_dp)
            endif
        enddo
    enddo

    call A%set_block(M,1,1)

    nullify(g,M)
    deallocate(nbrs)


    ! Create the graph and matrix for the triangles
    allocate(ll_graph::g)
    call g%init(n_triangles,n_triangles)

    do i=1,n_triangles
        call g%add_edge(i,i)
    enddo

    do i=1,n_triangles
        do k=1,3
            j = neigh(k,i)
            if (j/=-1) then
                call g%add_edge(i,j)
                call g%add_edge(j,i)
            endif
        enddo
    enddo

    call convert(g,'ell')

    allocate(ellpack_matrix::M)
    call M%init(n_triangles,n_triangles,'row',g)

    allocate(nbrs(g%max_degree))
    do i=1,n_triangles
        call g%neighbors(i,nbrs)
        do k=1,g%max_degree
            j = nbrs(k)
            if (j/=0 .and. j/=i) then
                call M%set_value(i,j,-1.0_dp)
                call M%add_value(i,i,1.0_dp)
            endif
        enddo
    enddo

    call A%set_block(M,2,2)

    nullify(g,M)


    ! Create the graph and matrix for the triangle-node connections
    allocate(ll_graph::g)
    call g%init(n_triangles,n_nodes)

    do n=1,n_triangles
        do k=1,3
            i = triangles(k,n)
            call g%add_edge(n,i)
        enddo
    enddo

    call convert(g,'ell')

    allocate(ellpack_matrix::M)
    call M%init(n_triangles,n_nodes,'row',g)

    call M%zero()

    call A%set_block(M,2,1)


    ! Create the graph and matrix for node-triangle interactions
    nullify(M)
    allocate(ellpack_matrix::M)
    call M%init(n_nodes,n_triangles,'col',g)

    call M%zero()

    call A%set_block(M,1,2)


    !----------------------------------------------------------------------!
    ! Test matrix-vector multiplication with the block matrix              !
    !----------------------------------------------------------------------!
    call x%init([n_nodes,n_triangles])
    call y%init([n_nodes,n_triangles])
    y%val = 0.0_dp

    do i=1,n_nodes
        call x%set_value([i,1],1.0_dp)
    enddo

    do n=1,n_triangles
        call x%set_value([i,2],1.0_dp)
    enddo

    call A%matmul(x,y)


end program block_matrix_tests

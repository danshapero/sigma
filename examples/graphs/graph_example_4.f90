!==========================================================================!
program graph_example_4                                                    !
!==========================================================================!
!==== Example program demonstrating graph edge iterators.              ====!
!====   The use of edge iterators is demonstrated by constructing a    ====!
!==== random graph using the Watts-Strogatz model for small-world      ====!
!==== networks.                                                        ====!
!====   Graph edge iterators provide a way of sequentially accessing   ====!
!==== all the edges of a graph with a well-defined interface, so that  ====!
!==== the user of a graph need not know precisely how that graph is    ====!
!==== implemented. While this operation can seem confusing or unusual  ====!
!==== to the un-initiated user, it has proved *indispensible* for the  ====!
!==== design of a robust sparse matrix class.                          ====!
!==========================================================================!

use sigma

implicit none

    ! graphs
    type(ll_graph) :: g, g_ring

    ! probability and random numbers
    real(dp) :: p, z, w

    ! integer indices
    integer :: i, j, k, nn

    ! variables for getting the neighbors of a vertex
    integer :: d
    integer, allocatable :: neighbors(:)

    ! variables for iterating through the edges of a graph
    integer :: n, num_batches, num_returned, edges(2,batch_size)
    type(graph_edge_cursor) :: cursor

    ! variables for computing the clustering coefficient of a graph
    integer :: k1, k2
    integer, allocatable :: adj(:,:)
    real(dp) :: local_clustering, clustering


    nn = 2048
    d = ceiling( log(1.0_dp * nn) / log(2.0) )

    call init_seed()


    !----------------------------------------------------------------------!
    ! Make the graph of a regular ring lattice                             !
    !----------------------------------------------------------------------!
    call g_ring%init(nn)

    do i = 1, nn
        do k = 1, d / 2
            j = mod(i + k - 1, nn) + 1
            call g_ring%add_edge(i, j)
            call g_ring%add_edge(j, i)
        enddo
    enddo

    write(*,10) nn, d
10  format('Generated a regular ring graph with 'i5,' vertices, degree ',i3,'.')



    !----------------------------------------------------------------------!
    ! Randomly rewire the edges of the ring lattice to generate a          !
    ! Watts-Strogatz graph                                                 !
    !----------------------------------------------------------------------!
    call g%init(nn)

    p = 0.25

    cursor = g_ring%make_cursor()

    num_batches = (cursor%final - cursor%start) / batch_size + 1
    do n = 1, num_batches
        ! Get a chunk of edges from the ring lattice
        call g_ring%get_edges(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            ! For each edge,
            i = edges(1, k)
            j = edges(2, k)

            if (j > i) then
                ! Flip a coin
                call random_number(z)

                ! If the coin turn up heads,
                if (z < p) then
                    ! Randomly select another vertex j, provided that j
                    ! is not equal to i and is not already connected to i
                    ! in the graph g
                    call random_number(w)
                    j = int(w * nn) + 1
                    do while(j == i .or. j > nn .or. g%connected(i, j))
                        call random_number(w)
                        j = int(w * nn) + 1
                    enddo
                endif

                ! Whether or not j has been switched to another vertex,
                ! connect i,j in g.
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo
    enddo

    write(*,20) p
20  format('Rewired all edges of ring graph with probability ', f8.6, '.')



    !----------------------------------------------------------------------!
    ! Compute the clustering coefficient of the Watts-Strogatz graph       !
    !----------------------------------------------------------------------!
	d = g%max_degree()
    allocate(neighbors(d), adj(d, d))

    clustering = 0.0_dp

    ! For each vertex k,
    do k = 1, nn
        ! find all the neighbors of k
        d = g%degree(k)
        call g%get_neighbors(neighbors, k)

        ! Fill the adjacency matrix of the sub-graph of g consisting of
        ! vertex k and all its neighbors
        adj = 0
        local_clustering = 0.0

        do k2 = 1, d
            j = neighbors(k2)
            do k1 = 1, d
                i = neighbors(k1)
                if (g%connected(i, j)) adj(k1, k2) = 1
            enddo
        enddo

        ! The local clustering coefficient is the number of entries in
        ! the local adjacency matrix, divided by the maximum possible
        ! number of entries, namely degree^2-degree.
        if (d > 1) local_clustering = sum(adj) / float(d * (d - 1))

        ! Add the local contribution to the global clustering coefficient
        clustering = clustering + local_clustering
    enddo

    ! The global clustering coefficient is the average of all the locals
    clustering = clustering / nn

    d = g_ring%max_degree()
    write(*,100) clustering, 3.0 * (d - 2) / (4.0 * (d - 1)), 1.0 * d / nn
100 format('Clustering for WS, ring & ER graphs: ',f8.6,', ',f8.6,', ',f8.6)


end program graph_example_4

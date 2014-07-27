!==========================================================================!
program graph_example_1                                                    !
!==========================================================================!
!==== Example program demonstrating the following graph operations:    ====!
!====       o  initialize a graph with a given number of vertices      ====!
!====       o  add edges into a graph                                  ====!
!====       o  check whether two vertices are connected                ====!
!====       o  find the degree of a vertex                             ====!
!==== These are illustrated by constructing a *random* graph using     ====!
!==== the *Erdos-Renyi* model: two vertices i, j are connected with    ====!
!==== some fixed probability p, independent of all other edges.        ====!
!==========================================================================!

use sigma

implicit none

    ! a linked-list graph and variables for generating it randomly
    type(ll_graph) :: g
    real(dp) :: z(512), p

    ! integer indices
    integer :: i, j, k, nn

    ! variables for estimating the probability that two vertices are
    ! connected
    integer :: num_connected

    ! variables for computing statistics on the degree of graph vertices
    integer :: d, min_degree, max_degree
    real(dp) :: avg_degree



    ! initialize a random seed and choose the probability of two vertices
    ! being connected
    call init_seed()
	nn = 512
    p = 6.25 / nn



    !----------------------------------------------------------------------!
    ! Set up a graph with edges chosen at random                           !
    !----------------------------------------------------------------------!

    ! Initialize g to be a graph with 512 vertices
    call g%init(nn)

    ! For each edge, pick some neighbors at random.
    do i = 1, nn
        ! Generate a whole mess of U[0,1] random numbers
        call random_number(z)

        do j = i + 1, nn
            ! For each j, if the random number z(j) drawn is less than p,
            ! make (i,j) connected in g.
            if (z(j) < p) then
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo
    enddo

    print *, 'Random graph generated with 512 vertices.'



    !----------------------------------------------------------------------!
    ! Estimate the probability that two nodes are connected                !
    !----------------------------------------------------------------------!
    num_connected = 0

    do k = 1, nn
        call random_number(z)

        i = int( nn*z(1)+1 )
        j = int( nn*z(2)+1 )

        if (g%connected(i, j)) then
            num_connected = num_connected + 1
        endif
    enddo

    print *, 'From a random sample of 1024 edges,'
    print *, num_connected, 'were connected.'



    !----------------------------------------------------------------------!
    ! Compute the minimum, maximum and average degree of the graph         !
    !----------------------------------------------------------------------!
    min_degree = nn
    max_degree = 0
    avg_degree = 0.0_dp

    do i = 1, nn
        d = g%degree(i)

        min_degree = min(d, min_degree)
        max_degree = max(d, max_degree)

        avg_degree = avg_degree + d
    enddo

    avg_degree = avg_degree / nn

    print *, 'Min/max/average degrees:', min_degree, max_degree, avg_degree



end program graph_example_1

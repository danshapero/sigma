!==========================================================================!
program matrix_example_1                                                   !
!==========================================================================!
!==== Example program demonstrating some basic matrix operations:      ====!
!====       o  initializing a matrix from a graph                      ====!
!====       o  set the entries of a matrix                             ====!
!====       o  get the value of a matrix entry                         ====!
!====       o  multiply a matrix by a vector                           ====!
!====       o  multiply the transpose of a matrix by a vector          ====!
!==== This program assumes the user has some familiarity with graph    ====!
!==== operations as a prerequisite.                                    ====!
!==== The matrix operations are illustrated by generating a random     ====!
!==== directed graph and simulating a random walk on this graph.       ====!
!==== The matrix generated is the transition probability matrix for    ====!
!==== underlying graph, assuming that the random walker is equally     ====!
!==== likely to go to any neighbor.                                    ====!
!==========================================================================!

use sigma

implicit none

    ! a linked-list graph and variables for generating it randomly
    class(graph_interface), pointer :: g
    real(dp) :: z(512), p, entropy

    ! a sparse matrix
    class(sparse_matrix_interface), pointer :: A

    ! some integer indices
    integer :: i, j, k, nn

    ! variables for making the transition matrix corresponding to a Markov
    ! chain on g
    integer :: d
    integer, allocatable :: neighbors(:)

    ! some vectors
    real(dp), allocatable :: x(:), y(:)


    nn = 512

    ! initialize a random seed
    call init_seed()
    p = 7.0 / nn

    allocate(x(nn), y(nn))



    !----------------------------------------------------------------------!
    ! Set up a random graph g                                              !
    !----------------------------------------------------------------------!
    call choose_graph_type(g, "list of lists")
    call g%init(nn)

    do i = 1, nn
        call random_number(z)

        do j = 1, nn
            if (z(j) < p) call g%add_edge(i, j)
        enddo
    enddo

    do i = 1, nn
        if (g%degree(i) == 0) call g%add_edge(i, i)
    enddo

    d = g%max_degree()
    write(*,10) d
10  format('Done generating random graph; max degree: ',i4)


    !----------------------------------------------------------------------!
    ! Create a matrix with the connectivity structure of g                 !
    !----------------------------------------------------------------------!
    ! Initialize a sparse matrix with 512 rows and columns, in row-major 
    ! ordering, with g representing its connectivity
    call choose_matrix_type(A, "compressed sparse row")

    allocate(neighbors(d))

    ! For each vertex i,
    do i = 1, nn
        ! compute the degree d of i
        d = g%degree(i)

        ! Find all the neighbors of i.
        call g%get_neighbors(neighbors, i)

        ! For each neighbor j,
        do k = 1, d
            j = neighbors(k)
            call A%set_value(i, j, 1.0_dp / d)
        enddo
    enddo



    !----------------------------------------------------------------------!
    ! Check to make sure that all the rows sum to 1                        !
    !----------------------------------------------------------------------!
    x = 1.0_dp
    y = 0.0_dp

    ! Compute the matrix-vector product,  y = A*x
    call A%matvec(x,y)

    write(*,20) minval(y) - 1.0_dp, maxval(y) - 1.0_dp
20  format('Rows sum to 1 to within error (', e12.6, ',', e12.6, ')')



    !----------------------------------------------------------------------!
    ! Starting from an initial state where the random walk is sure to be   !
    ! at vertex 1, look at the probability after several steps.            !
    !----------------------------------------------------------------------!
    x = 0.0
    x(1) = 1.0

    do k = 1, 4 * nn
        call A%matvec_t(x, y)
        x = y / sum(y)
    enddo

    entropy = 0.0
    do i = 1, nn
        if (x(i) /= 0) entropy = entropy - x(i) * log(x(i))
    enddo

    write(*,30) 100*entropy / log(1.0_dp * nn)
30  format('After 2048 iterations, entropy relative maximum is ',f9.6,'%.')

    call A%matvec_t(x, y)
    write(*,40)
40  format('Residual of invariant distribution as left eigenvector')
    write(*,50) maxval(dabs(x - y))
50  format('   of transition matrix: ', e12.6)


end program matrix_example_1

!==========================================================================!
program graph_example_1                                                    !
!==========================================================================!
!==== Example program demonstrating how to initialize a graph and fill ====!
!==== it with some random edges, then do a few operations on it.       ====!
!==========================================================================!

use fempack

implicit none

    ! a linked-list graph and variables for generating it randomly
    type(ll_graph) :: g
    real(dp) :: z(512), p

    ! integer indices
    integer :: i,j,k

    ! variables for computing the degree of the graph
    integer :: d, min_degree, max_degree
    real(dp) :: avg_degree

    ! variables for doing a breadth-first search through the graph
    integer, allocatable :: neighbors(:)
    integer :: component_size
    type(dynamic_array) :: stack
    logical :: found(512)



    ! initialize a random seed and choose the probability of two vertices
    ! being connected
    call init_seed()
    p = 6.25/512



    !----------------------------------------------------------------------!
    ! Set up a graph                                                       !
    !----------------------------------------------------------------------!

    ! Initialize g to be a graph with 512 vertices
    call g%init(512)

    ! For each edge, pick some neighbors at random.
    do i=1,512
        ! Generate a whole mess of U[0,1] random numbers
        call random_number(z)

        do j=i+1,512
            ! For each j, if the random number z(j) drawn is less than p,
            ! make (i,j) connected in g.
            if ( z(j)<p ) then
                call g%add_edge(i,j)
                call g%add_edge(j,i)
            endif
        enddo
    enddo



    !----------------------------------------------------------------------!
    ! Compute the minimum, maximum and average degree of the graph         !
    !----------------------------------------------------------------------!

    min_degree = 512
    max_degree = 0
    avg_degree = 0.0_dp

    do i=1,512
        d = g%degree(i)

        min_degree = min(d, min_degree)
        max_degree = max(d, max_degree)

        avg_degree = avg_degree+d
    enddo

    avg_degree = avg_degree/512

    print *, 'Random graph generated with 512 vertices and',g%ne,'edges.'
    print *, 'Min/max/average degrees:',min_degree,max_degree,avg_degree
    print *, ' '



    !----------------------------------------------------------------------!
    ! Depth-first search the graph                                         !
    !----------------------------------------------------------------------!

    allocate(neighbors(max_degree))
    found = .false.

    ! initialize a stack to store the edges that need visiting
    call stack%init()

    ! enqueue vertex #1
    call stack%push(1)
    found(1) = .true.

    ! continue until the queue is empty
    do while( stack%length>0 )
        ! pop a vertex i off the front of the queue
        i = stack%pop()

        ! find the degree of that vertex
        d = g%degree(i)

        ! find all the neighbors of vertex i
        call g%neighbors(i,neighbors)

        ! for each of those neighbors j,
        do k=1,d
            j = neighbors(k)
            ! check if vertex j has been visited yet.
            if (.not.found(j)) then
                ! if not, mark it as having been visited and enter it
                ! into the queue
                found(j) = .true.
                call stack%push(j)
            endif
        enddo
    enddo



    !----------------------------------------------------------------------!
    ! Find the size of the connected component of vertex 1                 !
    !----------------------------------------------------------------------!

    component_size = 0
    do i=1,512
        if (found(i)) component_size = component_size+1
    enddo

    print *, 'Depth-first search of graph complete.'
    if (component_size==512) then
        print *, 'Graph is connected! Hooray!'
    else
        print *, 'Graph is not connected! Oh bother!'
        print *, 'Size of connected component of vertex #1:',component_size
    endif


end program graph_example_1

!==========================================================================!
program graph_example_3                                                    !
!==========================================================================!
!==== Example program demonstrating more graph operations:             ====!
!====       o  deleting edges                                          ====!
!====       o  compressing graph storage                               ====!
!==== These are illustrated by way of *bond percolation*, a classic    ====!
!==== problem of physics. One considers an integer lattice in d-       ====!
!==== dimensional space, and every edge of the lattice is removed      ====!
!==== with probability p. What is the probability that there is a      ====!
!==== path through the lattice?                                        ====!
!==========================================================================!

use sigma

implicit none

    ! a CS-graph representing the lattice
    type(cs_graph) :: g
    real(dp) :: z(2), p

    ! integer indices
    integer :: i,j,k,x,y,xs,ys,nx,ny

    ! variables for doing a depth-first search through the graph
    integer :: components, d
    integer, allocatable :: neighbors(:), component(:,:)
    type(dynamic_array) :: stack

    ! variables for checking whether the graph succeeds at percolating
    integer :: xt, xb, comp, num_percolated
    integer, allocatable :: top_components(:), bottom_components(:)


    nx = 64
    ny = 64

    ! initialize a random seed and choose the probability of two vertices
    ! being connected
    call init_seed()
    p = 0.5



    !----------------------------------------------------------------------!
    ! Set up a graph representing the integer lattice on a grid of         !
    ! nx x ny points.                                                      !
    !----------------------------------------------------------------------!
    call g%init(nx*ny,degree=4)

    ! Make every (x,y) connected to (x,y+1), (x,y-1), (x+1,y) & (x-1,y)
    do x = 1,nx
        do y = 1,ny
            i = ny*(x-1)+y

            ! Add edge from (x,y) to (x,y+1) and vice versa
            j = ny*(x-1)+mod(y,ny)+1
            call g%add_edge(i,j)
            call g%add_edge(j,i)

            ! (x,y) to (x+1,y) and v.v.
            j = ny*mod(x,nx)+y
            call g%add_edge(i,j)
            call g%add_edge(j,i)
        enddo
    enddo

    ! Delete all the connections between (x,1) and (x,1024)
    do x=1,nx
        i = ny*(x-1)+1
        j = ny*(x-1)+ny

        call g%delete_edge(i,j)
        call g%delete_edge(j,i)
    enddo

    ! Delete all the connections between (1,y) and (1024,y)
    do y=1,ny
        i = y
        j = ny*(nx-1)+y

        call g%delete_edge(i,j)
        call g%delete_edge(j,i)
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
            i = ny*(x-1)+y

            ! If the first number is less than p, delete the connection
            ! to (x,y+1)
            j = ny*(x-1)+mod(y,ny)+1
            if (z(1)<p) then
                call g%delete_edge(i,j)
                call g%delete_edge(j,i)
            endif

            ! If the second number is less than p, delete the connection
            ! to (x+1,y)
            j = ny*mod(x,nx)+y
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
    ! Try to determine if the top of the lattice is still connected to     !
    ! the bottom by performing a whole mess of graph traversals.           !
    !----------------------------------------------------------------------!
    ! First, assume that the number of components is 0
    components = 0

    ! Make an array `components` which will indicate which connected
    ! component of the graph each node belongs to
    allocate(component(nx,ny),neighbors(g%max_degree))
    component = 0

    do while(minval(component)==0)
        call stack%init()

        ! Increment the number of components
        components = components+1

        ! Find a node which has not been assigned a component yet
        do y=1,ny
            do x=1,nx
                if (component(x,y)==0) then
                    xs = x
                    ys = y
                endif
            enddo
        enddo

        ! Compute the linear index i corresponding to the point (x,y)
        i = nx*(xs-1)+ys

        ! Indicate that i belongs to the current component
        component(xs,ys) = components

        ! Push i onto the stack
        call stack%push(i)

        ! Depth-first search g starting from i
        do while(stack%length>0)
            ! pop a vertex off the stack
            i = stack%pop()

            ! find the degree of that vertex
            d = g%degree(i)

            ! find all the neighbors of vertex i
            call g%neighbors(i,neighbors)

            ! for each of those neighbors j,
            do k=1,d
                j = neighbors(k)

                ! find the (x,y) coordinate of j
                y = mod(j,nx)+1
                x = (j-y)/nx+1

                ! check if vertex j has been visited yet
                if (component(x,y)==0) then
                    ! if not, mark it as belonging to the current component
                    component(x,y) = components

                    ! then push it onto the stack
                    call stack%push(j)
                endif
            enddo
        enddo

        ! clear the stack
        call stack%free()
    enddo

    write(*,30) components
30  format('Done traversing graph;       ',i8,' connected components.')



    !----------------------------------------------------------------------!
    ! Determine whether there is a connected component containing a point  !
    ! from both the top and the bottom of the lattice.                     !
    !----------------------------------------------------------------------!
    allocate(top_components(nx), bottom_components(nx))

    top_components = component(:,1)
    bottom_components = component(:,ny)

    num_percolated = 0
    do xt=1,nx
        comp = top_components(xt)

        do xb=1,nx
            if (bottom_components(xb)==comp) then
                num_percolated = num_percolated+1
                exit
            endif
        enddo
    enddo

    write(*,40) (100.0_dp*num_percolated)/nx
40  format(F12.6,'% of top sites percolated to the bottom.')


end program graph_example_3

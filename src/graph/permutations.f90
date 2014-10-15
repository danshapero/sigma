!==========================================================================!
!==========================================================================!
module permutations                                                        !
!==========================================================================!
!==========================================================================!
!==== This module contains subroutines for re-ordering graphs to       ====!
!==== accelerate sparse matrix computations, e.g.the Cuthill-McKee     ====!
!==== ordering, multicolor orderings, independent set orderings, etc.  ====!
!==========================================================================!
!==========================================================================!

use graph_interfaces
use types, only: dynamic_array, circular_array

implicit none

contains



!--------------------------------------------------------------------------!
subroutine breadth_first_search(g, p)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph_interface), intent(in) :: g
    integer, intent(out) :: p(:)
    ! local variables
    integer :: i, j, k, d, num
    integer, allocatable :: neighbors(:)
    type(circular_array) :: queue

    ! Allocate the neighbors array
    d = g%max_degree()
    allocate(neighbors(d))

    num = 0

    ! Initialize a queue, and enqueue vertex #1
    call queue%init(capacity = 64, min_capacity = 16)
    call queue%enqueue(1)

    ! Mark all nodes as unvisited
    p = -1

    ! So long as the queue is not empty,
    do while(queue%length > 0)
        num = num + 1

        ! Pop the first value i off
        i = queue%pop()

        ! Label i according to the order in which it was visited
        p(i) = num

        ! Find the degree and all neighbors of node `i`
        d = g%degree(i)
        call g%get_neighbors(neighbors, i)

        ! For each neighbor,
        do k = 1, d
            j = neighbors(k)
            
            ! if that neighbor has not been visited,
            if (p(j)==-1) then
                ! put it into the queue
                call queue%enqueue(j)

                ! and label it 0 to indicate that it has been queued.
                p(j) = 0
            endif
        enddo   ! End of loop over `k`

    enddo   ! End of while loop

    call queue%free()
    deallocate(neighbors)

end subroutine breadth_first_search



!--------------------------------------------------------------------------!
subroutine greedy_coloring(g, colors)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph_interface), intent(in) :: g
    integer, intent(out) :: colors(:)
    ! local variables
    integer :: i, j, k, d, used_colors, color, min_occupancy
    integer, allocatable , dimension(:) :: neighbors, neighbor_colors, &
                                            & color_totals
    type(circular_array) :: queue

    d = g%max_degree()
    allocate(neighbors(d), neighbor_colors(d + 1), color_totals(d + 1))

    ! Initialize a queue
    call queue%init(capacity = 32, min_capacity = 16)
    call queue%enqueue(1)

    ! Set all colors to a negative value to indicate they have not been
    ! visited or enqueued yet, and set color(1) to 0 to show that it's been
    ! queued up
    colors    = -1
    colors(1) = 0

    ! Set the tally of how many nodes there are of each color
    color_totals = 0
    used_colors  = 0

    ! Keep iterating until the queue is empty
    do while (queue%length > 0)
        neighbor_colors = 0

        ! Pop the first value i off the queue
        i = queue%pop()

        ! Find all the neighbors of node i
        d = g%degree(i)
        call g%get_neighbors(neighbors, i)

        ! Tally up how many neighbors there are of each color
        do k = 1, d
            j = neighbors(k)
            color = colors(j)
            if (color > 0) then
                neighbor_colors(color) = neighbor_colors(color) + 1
            elseif (color == -1) then
                call queue%enqueue(j)
                colors(j) = 0
            endif
        enddo

        color = 0
        min_occupancy = g%n + 1
        ! Find which neighboring color has the fewest nodes assigned to it
        do k = 1, used_colors
            if (color_totals(k) > 0 &
                    & .and. color_totals(k) < min_occupancy &
                    & .and. neighbor_colors(k) == 0) then
                color = k
                min_occupancy = color_totals(k)
            endif
        enddo

        ! If the color chosen is still zero, that means that, of all the
        ! colors we're already using, there is some neighbor that already
        ! has that color. This meas that node i needs to be given a brand
        ! new color we haven't used before.
        if (color == 0) then
            color = used_colors + 1
            used_colors = used_colors + 1
        endif

        colors(i) = color
        color_totals(color) = color_totals(color) + 1
    enddo

    call queue%free()
    deallocate(neighbors, neighbor_colors, color_totals)

end subroutine greedy_coloring



!--------------------------------------------------------------------------!
subroutine greedy_color_ordering(g, p, ptrs, num_colors)                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph_interface), intent(in) :: g
    integer, intent(out) :: p(:), ptrs(:), num_colors
    ! local variables
    integer :: i, k, d, color
    integer, allocatable :: added(:)

    ! Allocate an array to keep track of how many vertices of each color
    ! have been added to the permutation
    d = g%max_degree()
    allocate(added(d + 1))

    ! Assign colors to all the nodes and put them in the array p
    call greedy_coloring(g, p)

    ! Find the total number of colors used
    num_colors = maxval(p)

    ! Make an array for the starting point in the permutation of each color
    ptrs = 0
    do i = 1, g%n
        ptrs(p(i) + 1) = ptrs(p(i) + 1) + 1
    enddo

    ptrs(1) = 1
    do color = 1, num_colors
        ptrs(color + 1) = ptrs(color + 1) + ptrs(color)
    enddo

    ! Build the permutation
    added = 0
    do i = 1, g%n
        color = p(i)
        k = added(color)
        p(i) = ptrs(color) + k

        added(color) = added(color) + 1
    enddo

end subroutine greedy_color_ordering




end module permutations

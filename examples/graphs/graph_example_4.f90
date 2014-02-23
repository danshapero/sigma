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

    type(ll_graph) :: g
    real(dp) :: p
    real(dp), allocatable :: z(:)



end program graph_example_4

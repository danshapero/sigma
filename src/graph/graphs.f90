module graphs

use graph_interface
use ll_graphs
use coo_graphs
use cs_graphs
use ellpack_graphs

implicit none


! Parameter giving the number of different graph formats available.
! Test programs iterate over every graph type to check correctness of each
! graph implementation.
integer, parameter :: num_graph_types = 4



interface choose_graph_type
    module procedure choose_graph_type_by_int, choose_graph_type_by_name
end interface

interface convert_graph_type
    module procedure convert_graph_type_by_int, convert_graph_type_by_name
end interface

contains




!--------------------------------------------------------------------------!
subroutine choose_graph_type_by_name(g, frmt)                              !
!--------------------------------------------------------------------------!
!     Take in a polymorphic graph pointer and allocate to a graph type     !
! according to a string specifying the name of the desried storage format. !
!--------------------------------------------------------------------------!
    class(graph), pointer, intent(inout) :: g
    character(len=*), intent(in) :: frmt

    nullify(g)

    select case( trim(frmt) )
        case("ll", "lol", "list oflists")
            allocate(ll_graph::g)

        case("coo","coordinate")
            allocate(coo_graph::g)

        case("cs","compressed sparse","harwell-boeing","csr","csc")
            allocate(cs_graph::g)

        case("ellpack","itpack")
            allocate(ellpack_graph::g)
    end select

end subroutine choose_graph_type_by_name



!--------------------------------------------------------------------------!
subroutine choose_graph_type_by_int(g, t)                                  !
!--------------------------------------------------------------------------!
!     Take in a polymorphic graph pointer and allocate it to a graph type  !
! according to an integer input:                                           !
!     1: list-of-lists graph                                               !
!     2: coordinate graph                                                  !
!     3: compressed sparse graph                                           !
!     4: ellpack/itpack graph                                              !
! This subroutine is mostly to facilitate testing, so that we can easily   !
! iterate through every graph type.                                        !
!--------------------------------------------------------------------------!
    class(graph), pointer, intent(inout) :: g
    integer, intent(in) :: t

    nullify(g)

    select case(t)
        case(1)
            allocate(ll_graph::g)

        case(2)
            allocate(coo_graph::g)

        case(3)
            allocate(cs_graph::g)

        case(4)
            allocate(ellpack_graph::g)
    end select

end subroutine choose_graph_type_by_int



!--------------------------------------------------------------------------!
subroutine convert_graph_type_by_name(g, frmt)                             !
!--------------------------------------------------------------------------!
    ! input/output variables 
    class(graph), pointer, intent(inout) :: g
    character(len=*), intent(in) :: frmt
    ! local variables
    class(graph), pointer :: h

    h => g
    nullify(g)

    call choose_graph_type(g, frmt)
    call g%copy(h)

    call h%destroy()
    deallocate(h)

end subroutine convert_graph_type_by_name



!--------------------------------------------------------------------------!
subroutine convert_graph_type_by_int(g, t)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), pointer, intent(inout) :: g
    integer, intent(in) :: t
    ! local variables
    class(graph), pointer :: h

    h => g
    nullify(g)

    call choose_graph_type(g, t)
    call g%copy(h)

    call h%destroy()
    deallocate(h)

end subroutine convert_graph_type_by_int




end module graphs


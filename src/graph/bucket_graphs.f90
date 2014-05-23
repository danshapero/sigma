module bucket_graphs

use graph_interface
use types, only: dynamic_array

implicit none


!--------------------------------------------------------------------------!
type, extends(graph) :: bucket_graph                                       !
!--------------------------------------------------------------------------!
    integer :: num_buckets, bucket_size
    integer, allocatable :: edges(:,:), owner(:)
    type(dynamic_array), allocatable :: buckets(:)
    type(dynamic_array) :: freed_bucket_stack
contains
    !--------------
    ! Constructors
    procedure :: init_const_degree => bucket_init_const_degree
    procedure :: init_variable_degree => bucket_init_variable_degree
    procedure :: copy => bucket_copy

    !-----------
    ! Accessors
    procedure :: degree => bucket_degree
    procedure :: get_neighbors => bucket_get_neighbors
    procedure :: connected => bucket_connected
    procedure :: find_edge => bucket_find_edge

    !---------------
    ! Edge iterator
    procedure :: make_cursor => bucket_make_cursor
    procedure :: get_edges => bucket_get_edges

    !----------
    ! Mutators
    procedure :: add_edge => bucket_add_edge
    procedure :: delete_edge => bucket_delete_edge
    procedure :: left_permute => bucket_left_permute
    procedure :: right_permute => bucket_right_permute
    procedure :: compress => bucket_compress
    procedure :: decompress => bucket_decompress

    !-------------
    ! Destructors
    procedure :: destroy => bucket_destroy

    !--------------------------
    ! Testing, debugging & I/O
    procedure :: dump_edges => bucket_dump_edges
end type bucket_graph






contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine bucket_init_const_degree(g,n,m,degree)                          !
!--------------------------------------------------------------------------!
    class(bucket_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, degree


end subroutine bucket_init_const_degree



!--------------------------------------------------------------------------!
subroutine bucket_init_variable_degree(g,n,m,degrees)                      !
!--------------------------------------------------------------------------!
    class(bucket_graph), intent(inout) :: g
    integer, intent(in) :: n, degrees(:)
    integer, intent(in), optioanl :: m


end subroutine bucket_init_variable_degree



!--------------------------------------------------------------------------!
subroutine bucket_init_copy(g,h)                                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bucket_graph), intent(inout) :: g
    class(graph), intent(in)           :: h
    ! local variables
    integer :: i,j,k,n,num_returned,num_blocks,edges(2,64)
    type(graph_edge_cursor) :: cursor


end subroutine bucket_init_copy




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function bucket_degree(g,i) result(d)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bucket_graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: d
    ! local variables
    integer :: j,k,l,bucket

    d = 0
    do k=1,g%buckets(i)%length
        bucket = g%buckets(i)%get_entry(k)
        do l=1,g%bucket_size
            j = g%edges(l,bucket)
            if (j/=0) d = d+1
        enddo
    enddo

end function bucket_degree



!--------------------------------------------------------------------------!
subroutine bucket_get_neighbors(g,neighbors,i)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bucket_graph), intent(inout) :: g
    integer, intent(out) :: neighbors(:)
    integer, intent(in) :: i
    ! local variables
    integer :: j,k,l,bucket,next

    neighbors = 0
    next = 0
    do k=1,g%buckets(i)%length
        bucket = g%buckets(i)%get_entry(k)
        do l=1,g%bucket_size
            j = g%edges(l,bucket)
            if (j/=0) then
                next = next+1
                neighbors(next) = j
            endif
        enddo
    enddo

end subroutine bucket_get_neighbors



!--------------------------------------------------------------------------!
function bucket_connected(g,i,j)                                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bucket_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: bucket_connected
    ! local variables
    integer :: k,l

    bucket_connected = .false.

    do k=1,g%buckets(i)%length
        bucket = g%buckets(i)%get_entry(k)
        do l=1,g%bucket_size
            j = g%edges(l,bucket)
            if (j==i) bucket_connected = .true.
        enddo
    enddo

end function bucket_connected



!--------------------------------------------------------------------------!
function bucket_find_edge(g,i,j)                                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bucket_graph), intent(in) :: g
    integer, intent(in) :: i,j
    integer :: bucket_find_edge
    ! local variables
    integer :: k,l

    bucket_find_edge = -1

    do k=1,g%buckets(i)%length
        bucket = g%buckets(i)%get_entry(k)
        do l=1,g%bucket_size
            j = g%edges(l,bucket)
            if (j==i) k = g%bucket_size*(bucket-1)+l
        enddo
    enddo

end function bucket_find_edge




!==========================================================================!
!==== Edge iterator                                                    ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function bucket_make_cursor(g,thread) result(cursor)                       !
!--------------------------------------------------------------------------!
    class(bucket_graph), intent(in) :: g
    integer, intent(in) :: thread
    type(graph_edge_cursor) :: cursor

    cursor%start = 1
    cursor%final = g%bucket_size*g%num_buckets

end function bucket_make_cursor



!--------------------------------------------------------------------------!
subroutine bucket_get_edges(g,edges,cursor,num_edges,num_returned)         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bucket_graph), intent(inout) :: g
    integer, intent(out) :: edges(2,num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    ! local variables
    integer :: i,k,l,num_added,num_from_this_row

    ! Set up the returned edges to be 0
    edges = 0

    ! Count how many edges we're actually going to return
    num_returned = min(num_edges,cursor%final-cursor%current)

    ! 

end subroutine bucket_get_edges




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine bucket_add_edge(g,i,j)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bucket_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k,l,m,bucket
    logical :: added

    if (.not.g%mutable) then
        print *, 'Attempted to add an edge to immutable bucket graph'
        print *, 'Terminating'
        call exit(1)
    else
        if (.not.g%connected(i,j)) then
            added = .false.

            do k=1,g%buckets(i)%length
                bucket = g%buckets(i)%get_entry(k)
                do l=1,g%bucket_size
                    m = g%edges(l,bucket)
                    if (m==0) then
                        g%edges(l,bucket) = j
                        added = .true.
                        exit
                    endif
                enddo
            enddo

            if (.not.added) then

            endif
        endif
    endif

end subroutine bucket_add_edge



!--------------------------------------------------------------------------!
subroutine bucket_delete_edge(g,i,j)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bucket_graph), intent(inout) :: g
    integer :: i,j
    ! local variables
    integer :: k,l,bucket

    if (.not.g%mutable) then
        print *, 'Attempted to delete an edge from immutable bucket graph'
        print *, 'Terminating'
        call exit(1)
    else
        if (g%connected(i,j)) then

        endif
    endif

end subroutine bucket_delete_edge



end module bucket_graphs

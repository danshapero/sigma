module conversions

! Use all the graph modules
use graphs

implicit none

contains


!--------------------------------------------------------------------------!
subroutine convert(g,storage_format)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), pointer, intent(inout) :: g
    character(len=*), intent(in), optional :: storage_format
    ! local variables
    class(graph), pointer :: gc
    type(graph_edge_cursor) :: cursor
    integer :: num_nbrs(g%n),i,j,k
    integer :: n, num_batches, num_returned, edges(2,batch_size)
    character(len=3) :: str_fmt

    ! Ascertain the right storage format for g
    if (present(storage_format)) then
        str_fmt = storage_format
    else
        str_fmt = 'cs '
    endif

    ! Make a temporary pointer to g and nullify the original pointer to g
    gc => g
    nullify(g)

    ! Determine how many neighbors each node has
    num_nbrs = 0
    edges = 0
    cursor = gc%make_cursor()
    num_batches = (cursor%final-cursor%start+1)/batch_size+1

    do n=1,num_batches
        call gc%get_edges(edges,cursor,batch_size,num_returned)

        do k=1,num_returned
            i = edges(1,k)
            j = edges(2,k)

            if (i/=0 .and. j/=0) then
                num_nbrs(i) = num_nbrs(i)+1
            endif
        enddo
    enddo

    ! Allocate and initialize the new format for g
    select case(trim(str_fmt))
        case('ll')
            allocate(ll_graph::g)
        case('cs')
            allocate(cs_graph::g)
        case('ell')
            allocate(ellpack_graph::g)
        case('coo')
            allocate(coo_graph::g)
    end select

    call g%init(gc%n,gc%m,num_nbrs)

    ! Add in all the edges to the new format
    cursor = gc%make_cursor()

    do n=1,num_batches
        call gc%get_edges(edges,cursor,batch_size,num_returned)

        do k=1,num_returned
            i = edges(1,k)
            j = edges(2,k)

            if (i/=0 .and. j/=0) then
                call g%add_edge(i,j)
            endif
        enddo
    enddo

    ! Free up the space used for the original graph format
    call gc%destroy()
    deallocate(gc)

end subroutine convert



end module conversions

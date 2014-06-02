!--------------------------------------------------------------------------!
program graph_tests_2                                                      !
!--------------------------------------------------------------------------!
!     This program is for testing whether the graph copy constructor works !
! properly. To do this, a random linked-list graph `h` is generated. A     !
! graph `g` is copied from `h` for each graph type, and the two are        !
! checked for isomorphism.                                                 !
!--------------------------------------------------------------------------! 

use sigma

implicit none

    class(graph), pointer :: g, h
    integer :: i, j, k, l, n, degree, test1, test2, min_degree, ord(2)
    real(dp) :: z(16)
    logical :: trans
    ! variables for graph edge iteration
    integer :: num_blocks, num_returned, edges(2,batch_size)
    type(graph_edge_cursor) :: cursor
    ! command-line arguments
    character(len=16) :: arg
    logical :: verbose



    ! Get command line arguments to see if we're running in verbose mode
    verbose = .false.
    call getarg(1,arg)
    if (trim(arg)=="-v" .or. trim(arg)=="--verbose") then
        verbose = .true.
    endif

    ! Initialize a random seed
    call init_seed()


    !----------------------------------------------------------------------!
    ! Create a random sparse graph                                         !
    !----------------------------------------------------------------------!
    allocate(ll_graph::h)
    call h%init(512,512)

    do i=1,512
        ! Initialize the degree of node i to be 0
        degree = 0

        ! Make some random numbers
        call random_number(z)

        ! degree = number of samples greater than 1/2
        do k=1,16
            if (z(k)>0.5_dp) degree = degree+1
        enddo

        ! Make some more random numbers
        call random_number(z)

        do k=1,degree
            ! Get a random number j between 1 and 512
            j = min(int(z(k)*512)+1,512)

            ! Make (i,j) connected in h
            call h%add_edge(i,j)
        enddo
    enddo

    min_degree = h%max_degree

    select type(h)
        type is(ll_graph)
            do i=1,512
                min_degree = min( h%lists(i)%length, min_degree )
            enddo
    end select

    if (verbose) then
        print *, 'Random graph generated'
        print *, 'Min/max degree: ',min_degree,h%max_degree
        print *, ' '
    endif


    !----------------------------------------------------------------------!
    ! Copy the reference graph h to graphs in other formats                !
    !----------------------------------------------------------------------!
    if (verbose) print *, 'Testing copy constructor'

    do test1=1,4
        do test2=1,2
            ! Iterate over every graph format
            call choose_graph_type(g,test1)
            if (verbose) print *, '   Test', test1, test2

            ! Iterate over each graph orientation, i.e. test copying both h
            ! and h with all directed edges reversed
            trans = .false.
            ord = [1, 2]
            if (test2==2) then
                trans = .true.
                ord = [2, 1]
            endif

            call g%copy(h, trans)

            ! Iterate through all the edges of g and make sure they all are
            ! connected in h
            cursor = g%make_cursor(0)
            num_blocks = (cursor%final-cursor%start)/batch_size+1

            do n=1,num_blocks
                call g%get_edges(edges,cursor,batch_size,num_returned)

                do k=1,num_returned
                    i = edges(ord(1),k)
                    j = edges(ord(2),k)

                    if (i/=0 .and. j/=0) then
                        l = h%find_edge(i,j)
                        if (l==-1) then
                            print *, '    On test',test1,test2
                            print *, '    Vertices',i,j
                            print *, '    not connected in reference graph,'
                            print *, '    but are connected in copy graph.'
                            call exit(1)
                        endif
                    endif
                enddo
            enddo

            ! Iterate through all the edges of h and make sure they all are
            ! connected in g
            cursor = h%make_cursor(0)
            num_blocks = (cursor%final-cursor%start)/batch_size+1

            do n=1,num_blocks
                call h%get_edges(edges,cursor,batch_size,num_returned)

                do k=1,num_returned
                    i = edges(ord(1),k)
                    j = edges(ord(2),k)

                    if (i/=0 .and. j/=0) then
                        l = g%find_edge(i,j)
                        if (l==-1) then
                            print *, '    On test',test1,test2
                            print *, '    Vertices',i,j
                            print *, '    are connected in reference graph,'
                            print *, '    but not connected in copy graph.'
                            call exit(1)
                        endif
                    endif
                enddo
            enddo

            deallocate(g)
        enddo
    enddo


end program graph_tests_2

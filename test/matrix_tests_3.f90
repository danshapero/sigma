!--------------------------------------------------------------------------!
program matrix_tests_3                                                     !
!--------------------------------------------------------------------------!
!     This program tests setting matrix entries where there was no space   !
! pre-allocated for them before.                                           !
!--------------------------------------------------------------------------!

use sigma

implicit none

    type(sparse_matrix) :: A
    class(graph), pointer :: g, h
    integer :: test
    integer :: i,j,nn
    real(dp), allocatable :: x(:), y(:)
    ! command-line arguments
    character(len=16) :: arg
    logical verbose


    ! Get command line arguments to see if we're running in verbose mode
    verbose = .false.
    call getarg(1,arg)
    select case(trim(arg))
        case("-v")
            verbose = .true.
        case("-V")
            verbose = .true.
        case("--verbose")
            verbose = .true.
    end select


    nn = 64

    allocate(ll_graph::h)
    call h%init(nn,degree=3)
    do i=1,nn
        call h%add_edge(i,i)

        j = mod(i,nn)+1
        call h%add_edge(i,j)

        j = mod(i+nn-2,nn)+1
        call h%add_edge(i,j)
    enddo

    allocate(x(nn), y(nn))



    !----------------------------------------------------------------------!
    ! Construct graph for connectivity structure of sparse matrix          !
    !----------------------------------------------------------------------!
    do test=1,4
        if (verbose) print *, 'Test:',test

        call choose_graph_type(g,test)
        call g%init(h)


        !------------------------------------------------------------------!
        ! Make a matrix with g as its connectivity structure               !
        !------------------------------------------------------------------!
        call A%init(nn,nn,'row',g)
        call A%zero()

        do i=1,nn
            call A%set_value(i,i,5.0_dp)

            j = mod(i,nn)+1
            call A%set_value(i,j,-1.0_dp)
            call A%set_value(j,i,-1.0_dp)
        enddo

        if (A%g%max_degree/=3) then
            print *, 'Max degree of A%g should be 3 after initializing and'
            print *, 'filling circulant matrix.'
            print *, 'Value found:',A%g%max_degree
            print *, 'Terminating'
            call exit(1)
        endif

        if (A%g%connected(1,5) .or. A%g%connected(1,7)) then
            print *, 'A%g should not have 1,5 or 1,7 connected yet.'
            print *, 'Terminating.'
            call exit(1)
        endif

        x = 1.0_dp
        y = 0.0_dp
        call A%matvec(x,y)

        if ( minval(y)<3-1.0d-14 .or. maxval(y)>3+1.0d-14) then
            print *, 'All row sums of A should be 3.0;'
            print *, 'Range of row sums:',minval(y),maxval(y)
            print *, 'Terminating.'
            call exit(1)
        endif


        !------------------------------------------------------------------!
        ! Add some entries to A that weren't present in g                  !
        !------------------------------------------------------------------!
        do i=1,nn
            j = mod(i+3,nn)+1

            call A%set_value(i,j,-1.0_dp)
            call A%set_value(j,i,-1.0_dp)
        enddo

        if (A%g%max_degree/=5) then
            print *, 'Max degree of A%g should be 5 after adding 2 matrix '
            print *, 'entries to every row where they were not already pre-'
            print *, 'allocated. Value found:',A%g%max_degree
            print *, 'Terminating.'
            call exit(1)
        endif

        if (A%nnz/=5*nn) then
            print *, 'Number of non-zero entries of A after adding 2 matrix'
            print *, 'entries to every row where they were not already pre-'
            print *, 'allocated should be',5*nn
            print *, 'Value found:',A%nnz
            print *, 'Terminating.'
            call exit(1)
        endif

        x = 1.0_dp
        y = 0.0_dp

        call A%matvec(x,y)

        if (minval(y)<1-1.0d-14 .or. maxval(y)>1+1.0d-14) then
            print *, 'After adding new -1 entries, all row sums of A'
            print *, 'should be 1.0.'
            print *, 'Range of row sums:',minval(y),maxval(y)
            print *, 'Terminating.'
        endif


        ! Add some more entries
        do i=1,nn
            j = mod(i+5,nn)+1

            call A%set_value(i,j,-1.0_dp)
            call A%set_value(j,i,-1.0_dp)
        enddo

        ! Clear all data for the next test
        call A%destroy()
        call g%destroy()
        deallocate(g)
    enddo


end program matrix_tests_3

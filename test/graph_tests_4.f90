!--------------------------------------------------------------------------!
program graph_tests_4                                                      !
!--------------------------------------------------------------------------!
!    This program tests the graph union and product operations. These are  !
! used behind the scenese in constructing the sum and product of sparse    !
! matrices.                                                                !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! Matrices and graphs
    class(graph), pointer :: g, h1, h2, gr, hr
    integer, allocatable :: A(:,:), B1(:,:), B2(:,:)
    ! Assorted variables
    integer :: i, j, nn, test1, test2, test3, t1, t2
    logical :: tr1, tr2
    ! Random numbers
    real(dp) :: p
    real(dp), allocatable :: y(:), z(:)
    ! Command-line arguments
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


    ! Initialize a random seed
    call init_seed()
    nn = 64
    p = 8.0/nn


    !----------------------------------------------------------------------!
    ! Construct reference graphs from which all test graphs are copied     !
    !----------------------------------------------------------------------!
    allocate(ll_graph::gr)
    allocate(ll_graph::hr)

    call gr%init(nn)
    call hr%init(nn)

    allocate(y(nn),z(nn),A(nn,nn),B1(nn,nn),B2(nn,nn))
    do i=1,nn
        call random_number(z)
        call random_number(y)

        do j=1,nn
            if (z(j)<p) call gr%add_edge(i,j)
            if (y(j)<p) call hr%add_edge(i,j)
        enddo
    enddo


    !----------------------------------------------------------------------!
    ! Test graph union operation                                           !
    !----------------------------------------------------------------------!
    if (verbose) print *, '>> Testing graph union operation.'

    do test1=1,4
        if (verbose) print *, '    Test:',test1
        call choose_graph_type(h1,test1)
        call h1%init(gr)

        do test2=1,4
            if (verbose) print *, '        Test:',test2
            call choose_graph_type(h2,test2)
            call h2%init(hr)

            do test3=1,4
                if (verbose) print *, '            Test:',test3
                call choose_graph_type(g,test3)

                ! Check each possibility for transposing the input graphs
                ! h1, h2 of the graph union.
                do t1=1,2
                do t2=1,2
                    tr1 = (t1==1)
                    tr2 = (t2==1)

                    ! Compute the graph union
                    call graph_union(g, h1, h2, trans1=tr1, trans2=tr2)

                    ! Check that g is isomorphic to h1 + h2
                    call g%to_dense_graph(A)
                    call h1%to_dense_graph(B1, trans=tr1)
                    call h2%to_dense_graph(B2, trans=tr2)

                    A = A-(B1+B2)
                    if (maxval(A)/=0) then
                        print *, 'Graph union failed on test', &
                            & test1,test2,test3
                        print *, 'With transpositions', tr1, tr2
                        print *, 'Terminating.'
                        call exit(1)
                    endif

                    call g%destroy()
                enddo
                enddo

                deallocate(g)
            enddo

            call h2%destroy()
            deallocate(h2)
        enddo

        call h1%destroy()
        deallocate(h1)
    enddo



    !----------------------------------------------------------------------!
    ! Test graph product operation                                         !
    !----------------------------------------------------------------------!
    print *, ' '
    if (verbose) print *, '>> Testing graph product operation.'

    do test1=1,4
        if (verbose) print *, '    Test:',test1
        call choose_graph_type(h1,test1)
        call h1%init(gr)

        do test2=1,4
            if (verbose) print *, '        Test:',test2
            call choose_graph_type(h2,test2)
            call h2%init(hr)

            do test3=1,4
                if (verbose) print *, '            Test:',test3
                call choose_graph_type(g,test3)

                ! Check each possibility for transposing the input graphs
                ! h1, h2 of the graph union.
                do t1=1,2
                do t2=1,2
                    tr1 = (t1==1)
                    tr2 = (t2==1)

                    ! Compute the graph union
                    call graph_product(g, h1, h2, trans1=tr1, trans2=tr2)

                    ! Check that g is isomorphic to h1 + h2
                    call g%to_dense_graph(A)
                    call h1%to_dense_graph(B1, trans=tr1)
                    call h2%to_dense_graph(B2, trans=tr2)

                    A = A-matmul(B1,B2)
                    if (maxval(A)/=0) then
                        print *, 'Graph product failed on test', &
                            & test1,test2,test3
                        print *, 'With transpositions', tr1, tr2
                        print *, 'Terminating.'
                        call exit(1)
                    endif

                    call g%destroy()
                enddo
                enddo

                deallocate(g)
            enddo

            call h2%destroy()
            deallocate(h2)
        enddo

        call h1%destroy()
        deallocate(h1)
    enddo




end program graph_tests_4

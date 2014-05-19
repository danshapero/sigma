program solver_tests_2

use sigma

implicit none

    ! graph and matrix
    class(graph), pointer :: g
    type(sparse_matrix) :: A
    ! vectors
    real(dp), allocatable :: x(:), xc(:), b(:)
    ! direct solver
    class(cg_solver), pointer :: itsol
    class(sparse_ldu_solver), pointer :: ilu
    ! indices and miscellaneous
    integer :: i,j,k,d,nn
    integer, allocatable :: neighbors(:)
    real(dp) :: p, z
    ! command-line arguments
    character(len=16) :: arg
    logical :: verbose


    nn = 128
    p = log(1.0_dp*nn)/nn
    !call init_seed()

    !----------------------------------------------------------------------!
    ! Get command line arguments to see if we're running in verbose mode   !
    !----------------------------------------------------------------------!
    verbose = .false.
    call getarg(1,arg)
    if (trim(arg)=="-v" .or. trim(arg)=="--verbose") then
        verbose = .true.
    endif


    !----------------------------------------------------------------------!
    ! Create a random graph                                                !
    !----------------------------------------------------------------------!
    allocate(ll_graph::g)
    call g%init(nn)
    do i=1,nn
        call g%add_edge(i,i)

        do j=i+1,nn
            call random_number(z)

            if (z<p) then
                call g%add_edge(i,j)
                call g%add_edge(j,i)
            endif
        enddo
    enddo

    call g%compress()

    if (verbose) then
        print *, 'Done creating random connectivity graph G.'
        print *, 'Nodes:',nn
        print *, 'Edges:',g%ne/2
        print *, 'Maximum degree:',g%max_degree
    endif


    !----------------------------------------------------------------------!
    ! Make the Laplacian of the graph                                      !
    !----------------------------------------------------------------------!
    call A%init(nn,nn,'row',g)

    allocate(neighbors(g%max_degree))

    do i=1,nn
        call A%set_value(i,i,1.0_dp)

        call g%get_neighbors(neighbors,i)

        d = g%degree(i)
        do k=1,d
            j = neighbors(k)

            if (j/=i) then
                call A%set_value(i,j,-1.0_dp)
                call A%add_value(i,i,+1.0_dp)
            endif
        enddo
    enddo

    if (verbose) then
        print *, 'Done creating A = I + L(G).'
    endif


    !----------------------------------------------------------------------!
    ! Build a direct solver for the system                                 !
    !----------------------------------------------------------------------!
    allocate(ilu)
    call ilu%setup(A)

    allocate(itsol)
    call itsol%setup(A)


    !----------------------------------------------------------------------!
    ! Solve a linear system with the direct solver                         !
    !----------------------------------------------------------------------!
    allocate(x(nn), xc(nn), b(nn))
    x = 0.0_dp
    call random_number(xc)
    call A%matvec(xc,b)

    call ilu%solve(A,x,b)

    do i=1,nn
        if (isnan(x(i))) print *, i
    enddo

    

end program solver_tests_2

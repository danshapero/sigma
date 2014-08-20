!--------------------------------------------------------------------------!
program solver_test_incomplete_cholesky                                    !
!--------------------------------------------------------------------------!
!     This program tests using incomplete matrix factorizations, both as   !
! a stand-alone solver and as a preconditioner.                            !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph
    class(graph_interface), pointer :: g

    ! sparse matrix
    class(sparse_matrix), pointer :: A

    ! linear solver
    class(linear_solver), pointer :: solver, pc

    ! vectors
    real(dp), dimension(:), allocatable :: u, v, f, r, q

    ! neighbor list
    integer, allocatable :: nodes(:)

    ! integer indices
    integer :: i, j, k, n, d, nn

    ! random numbers
    real(dp) :: p, z

    ! error of linear solver
    real(dp) :: misfit

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose



    !----------------------------------------------------------------------!
    ! Get command line arguments to see if we're running in verbose mode   !
    !----------------------------------------------------------------------!
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



    !----------------------------------------------------------------------!
    ! Set the system size and initialize a random seed                     !
    !----------------------------------------------------------------------!

    nn = 128
    p = log(1.0_dp * nn) / log(2.0_dp) / nn

    call init_seed()



    !----------------------------------------------------------------------!
    ! Make a random graph `g`                                              !
    !----------------------------------------------------------------------!

    allocate( ll_graph :: g )
    call g%init(nn)

    do i = 1, nn
        call g%add_edge(i, i)

        do j = i + 1, nn
            call random_number(z)

            if (z < p) then
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo
    enddo

    d = g%max_degree()
    allocate(nodes(d))

    if (verbose) then
        print *, 'o Done generating random graph.'
        print *, '    Number of vertices:', nn
        print *, '    Number of edges:   ', g%ne
        print *, '    Max vertex degree: ', d
    endif



    !----------------------------------------------------------------------!
    ! Make a random matrix `A` on the graph `g`                            !
    !----------------------------------------------------------------------!

    A => sparse_matrix(nn, nn, g, 'row')
    call A%zero()

    do i = 1, nn
        call A%add_value(i, i, 1.0_dp)

        call g%get_neighbors(nodes, i)
        d = g%degree(i)
        do k = 1, d
            j = nodes(k)

            call random_number(z)
            if (j > i) then
                call A%set_value(i, j, -z)
                call A%add_value(i, i, +z)

                call A%set_value(j, i, -z)
                call A%add_value(j, j, +z)
            endif
        enddo
    enddo

    if (verbose) print *, 'o Done building random matrix.'



    !----------------------------------------------------------------------!
    ! Build incomplete Cholesky and CG solvers for `A`                     !
    !----------------------------------------------------------------------!

    solver => cg(1.d-16)
    pc => ldu(incomplete = .true., level = 0)

    call solver%setup(A)
    call pc%setup(A)

    if (verbose) print *, 'o Done building solver & preconditioner.'



    !----------------------------------------------------------------------!
    ! Create a random vector and generate a right-hand side from it        !
    !----------------------------------------------------------------------!

    allocate( u(nn), v(nn), f(nn), r(nn), q(nn) )

    u = 0.0_dp
    v = 0.0_dp
    f = 0.0_dp
    r = 0.0_dp
    q = 0.0_dp

    ! Make v random
    call random_number(v)

    ! Apply the the preconditioner once to smooth it over somewhat.
    ! A totally random vector is not very smooth and its dominant components
    ! will be in the span of `A`-eigenvectors corresponding to the way
    ! upper end of the spectrum.
    call A%matvec(v, q)
    r = v - q
    call pc%solve(A, v, r)

    ! Make the right-hand side vector `f` the product of `A` and `v`
    call A%matvec(v, f)

    u = 0.0_dp
    r = 0.0_dp



    !----------------------------------------------------------------------!
    ! Test using incomplete Cholesky as a solver                           !
    !----------------------------------------------------------------------!

    ! Since `u` is 0, the residual is initially equal to `f`
    r = f
    q = 0.0_dp

    do n = 1, 10 * nn
        ! Update `u`
        call pc%solve(A, q, r)
        u = u + q

        ! Update the residual `r`
        call A%matvec(u, q)
        r = f - q
    enddo

    misfit = maxval(dabs(u - v))

    if (misfit > 1.0e-14) then
        print *, 'Incomplete Cholesky failed to produce sufficiently'
        print *, 'accurate solution after', 10 * nn, 'iterations.'
        print *, 'Error:', misfit
        call exit(1)
    endif

    if (verbose) then
        print *, 'o Done solving with incomplete Cholesky.'
        print *, '    Error:', misfit
    endif



    !----------------------------------------------------------------------!
    ! Test using incomplete Cholesky as a preconditioner                   !
    !----------------------------------------------------------------------!

    u = 0.0_dp
    call solver%solve(A, u, f, pc)

    misfit = maxval(dabs(u - v))

    if (misfit > 1.0e-15) then
        print *, 'ILDU-preconditioned conjugate gradient method failed'
        print *, 'to produce sufficiently accurate solution.'
        print *, 'Error:', misfit
        call exit(1)
    endif

    if (verbose) then
        print *, 'o Done solving with ILDU-preconditioned CG.'
        print *, '    Error:', misfit
    endif





    call g%destroy()
    call A%destroy()
    call solver%destroy()
    call pc%destroy()
    deallocate(g, A, solver, pc)


end program solver_test_incomplete_cholesky

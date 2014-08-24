!--------------------------------------------------------------------------!
program solver_test_diffusion_1d                                           !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph
    class(graph_interface), pointer :: g

    ! sparse matrix
    type(ellpack_matrix) :: A

    ! linear solver
    class(linear_solver), pointer :: solver

    ! vectors
    real(dp), dimension(:), allocatable :: u, v, f

    ! integer indices
    integer :: i, nn

    ! mesh spacing
    real(dp) :: dx

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
    ! Create a matrix discretizing the operator - d^2 / dx^2               !
    !----------------------------------------------------------------------!
    nn = 127
    dx = 1.0_dp / (nn + 1)

    allocate(ll_graph :: g)
    call g%init(nn, nn)
    do i = 1, nn - 1
        call g%add_edge(i, i)
        call g%add_edge(i, i + 1)
        call g%add_edge(i + 1, i)
    enddo
    call g%add_edge(nn, nn)

    call convert_graph_type(g, "ellpack")

    call A%set_dimensions(nn, nn)
    call A%set_graph(g)
    call A%zero()

    do i = 1, nn - 1
        call A%set_value(i, i,     +2.0_dp)
        call A%set_value(i, i + 1, -1.0_dp)
        call A%set_value(i + 1, i, -1.0_dp)
    enddo
    call A%set_value(nn, nn, 2.0_dp)

    if (verbose) print *, 'Done creating matrix for 1d Laplace operator'



    !----------------------------------------------------------------------!
    ! Create a reference solution for the linear system                    !
    !----------------------------------------------------------------------!

    allocate( u(nn), v(nn), f(nn) )
    u = 0.0_dp
    f = 2.0 * dx**2

    do i = 1, nn
        v(i) = i * dx * (1.0_dp - i * dx)
    enddo

    if (verbose) print *, 'Done creating reference solution'



    !----------------------------------------------------------------------!
    ! Test that the solver gives the correct answer                        !
    !----------------------------------------------------------------------!

    solver => cg(1.d-16)
    call solver%setup(A)

    if (verbose) print *, 'Done setting up linear solver'

    call solver%solve(A, u, f)

    misfit = maxval(dabs(u - v))

    if (verbose) print *, 'Done solving system'

    if ( misfit > 1.0e-14 ) then
        print *, 'CG solver failed.'
        print *, 'Should have error <', 1.0e-14
        print *, 'Error found:', misfit
        call exit(1)
    endif

    if (verbose) print *, 'Error:', misfit



    call g%destroy()
    call A%destroy()
    call solver%destroy()

    deallocate(g, solver)

end program solver_test_diffusion_1d

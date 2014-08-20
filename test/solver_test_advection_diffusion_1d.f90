!--------------------------------------------------------------------------!
program solver_test_advection_diffusion_1d                                 !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph
    class(graph_interface), pointer :: g

    ! sparse matrix
    class(sparse_matrix), pointer :: A

    ! linear_solver
    class(linear_solver), pointer :: solver

    ! vectors
    real(dp), dimension(:), allocatable :: u, v, f

    ! integer indices
    integer :: i, nn

    ! position and mesh spacing
    real(dp) :: x, dx

    ! advection velocity
    real(dp) :: c

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
    ! Create a matrix discretizing the operator - d^2 / dx^2 + c * d / dx  !
    !----------------------------------------------------------------------!
    nn = 1024
    dx = 1.0_dp / (nn + 1)
    c = 0.5_dp

    allocate(ll_graph :: g)
    call g%init(nn, nn)
    do i = 1, nn - 1
        call g%add_edge(i, i)
        call g%add_edge(i, i + 1)
        call g%add_edge(i + 1, i)
    enddo
    call g%add_edge(nn, nn)

    A => sparse_matrix(nn, nn, g, 'row')

    do i = 1, nn - 1
        call A%set_value(i, i,     +2.0_dp)
        call A%set_value(i, i + 1, -1.0_dp + c * dx/2)
        call A%set_value(i + 1, i, -1.0_dp - c * dx/2)
    enddo
    call A%set_value(nn, nn, 2.0_dp)

    if (verbose) then
        print *, 'Done creating matrix for 1d advection/diffusion operator'
    endif



    !----------------------------------------------------------------------!
    ! Create a reference solution for the linear system                    !
    !----------------------------------------------------------------------!

    allocate( u(nn), v(nn), f(nn) )
    u = 0.0_dp
    f = 2.0_dp * dx**2

    do i = 1, nn
        x = i * dx
        v(i) = 2.0_dp * ( x - (exp(c * x) - 1) / (exp(c) - 1) )/c
    enddo

    if (verbose) print *, 'Done creating reference solution'



    !----------------------------------------------------------------------!
    ! Test that the linear solver gives the correct answer                 !
    !----------------------------------------------------------------------!

    solver => bicgstab(1.0d-12)
    call solver%setup(A)

    if (verbose) print *, 'Done setting up linear solver'

    call solver%solve(A, u, f)

    misfit = maxval(dabs(u - v))

    if (verbose) print *, 'Done solving system'

    if ( misfit > 1.0e-8 ) then
        print *, 'BiCG-Stab solver failed.'
        print *, 'Should have error <', 1.0e-8
        print *, 'Error found:', misfit
        call exit(1)
    endif

    if (verbose) print *, 'Error:', misfit



    call g%destroy()
    call A%destroy()
    call solver%destroy()

    deallocate(g, A, solver)



end program solver_test_advection_diffusion_1d

program poisson

use fempack
use fem
use omp_lib

implicit none

    ! variables for constructing the computational mesh
    character(len=64) :: filename
    class(graph), pointer :: g
    real(dp), allocatable :: x(:,:)
    integer, allocatable :: bnd(:), ele(:,:), mask(:)

    ! variables for the stiffness and mass matrices
    class(sparse_matrix), pointer :: A, B

    ! variables for solving the linear system
    class(iterative_solver), allocatable :: solver
    class(preconditioner), allocatable :: pc
    real(dp), allocatable :: u(:), f(:), r(:), z(:), p(:)

    ! other variables
    integer :: nn, ne, n, next
    ! integer, allocatable :: nbrs(:)


!--------------------------------------------------------------------------!
! Read in the triangular mesh and assemble the matrices                    !
!--------------------------------------------------------------------------!
    call get_environment_variable('FEMPACK',filename)
    filename = trim(filename)//'/examples/meshes/circle.1'
    call read_triangle_mesh(g,x,bnd,ele,filename)

    nn = g%n
    ne = g%ne

    allocate(mask(sum(bnd)))
    next = 0
    do n=1,nn
        if (bnd(n)==1) then
            next = next+1
            mask(next) = n
        endif
    enddo

    allocate(cs_matrix::A)
    allocate(cs_matrix::B)
    call A%init(nn,nn,'row',g)
    call B%init(nn,nn,'row',g)


!--------------------------------------------------------------------------!
! Fill in the stiffness and mass matrices                                  !
!--------------------------------------------------------------------------!
    call laplacian2d(A,x,ele)
    call mass2d(B,x,ele)

    call A%sub_matrix_add(B)

    do n=1,size(mask)
        call A%add_value(mask(n),mask(n),1.0d8)
    enddo


!--------------------------------------------------------------------------!
! Allocate vectors for the RHS & solution                                  !
!--------------------------------------------------------------------------!
    allocate(u(nn),f(nn),z(nn),r(nn),p(nn))

    f = 1.0_dp
    u = 0.0_dp
    call B%matvec(f,u)
    f = u
    u = 0.0_dp


!--------------------------------------------------------------------------!
! Solve the linear system                                                  !
!--------------------------------------------------------------------------!
    allocate(cg_solver::solver)
    allocate(jacobi_preconditioner::pc)

    call solver%init(nn)
    call pc%init(A)

    call solver%solve(A,u,f,pc)
    print *, 'CG iterations to solve system:',solver%iterations
    print *, 'Range of solution:',minval(u),maxval(u)



end program

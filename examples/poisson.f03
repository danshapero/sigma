program poisson

use fempack
use fem

implicit none


    ! variables for constructing the computational mesh
    class(graph), pointer :: g
    real(dp), allocatable :: x(:,:)
    integer, allocatable :: bnd(:), ele(:,:), mask(:)

    ! variables for the stiffness and mass matrices
    class(sparse_matrix), allocatable :: A, B

    ! variables for solving the linear system
    class(iterative_solver), allocatable :: solver
    class(preconditioner), allocatable :: pc
    real(dp), allocatable :: u(:), f(:), r(:), z(:), p(:)
    real(dp) :: alpha, beta, dpr, res2

    ! other variables
    integer :: nn, ne, n, next
    integer, allocatable :: nbrs(:)


!--------------------------------------------------------------------------!
! Read in the triangular mesh and assemble the matrices                    !
!--------------------------------------------------------------------------!
    call read_triangle_mesh(g,x,bnd,ele,'meshes/circle.1')

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

    allocate(csr_matrix::A)
    allocate(csr_matrix::B)
    call A%assemble(g)
    call B%assemble(g)


!--------------------------------------------------------------------------!
! Fill in the stiffness and mass matrices                                  !
!--------------------------------------------------------------------------!
    call laplacian2d(A,x,ele)
    call mass2d(B,x,ele)

    do n=1,size(mask)
        call A%add_value(mask(n),mask(n),1.0d8)
    enddo


!--------------------------------------------------------------------------!
! Allocate vectors for the RHS, solution and CG iteration                  !
!--------------------------------------------------------------------------!
    allocate(u(nn),f(nn),z(nn),r(nn),p(nn))

    f = 1.0_dp
    u = 0.0_dp
    call B%matvec(f,u)
    f = u
    !u = 0.0_dp

    allocate(cg_solver::solver)
    allocate(jacobi_preconditioner::pc)

    call solver%init(nn)
    call pc%init(A)

    call solver%solve(A,u,f,pc)
    print *, solver%iterations

    print *, minval(u),maxval(u)


!--------------------------------------------------------------------------!
! Write out the solution to a file                                         !
!--------------------------------------------------------------------------!
    open(unit=10,file='u.txt')
    do n=1,nn
        write(10,*) u(n)
    enddo
    close(10)

end program

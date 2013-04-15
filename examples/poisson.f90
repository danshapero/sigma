program poisson

    use mesh_mod
    use linear_algebra_mod
    use fem_mod
    use netcdf

    implicit none

    ! command line arguments
    character(len=32) :: meshname,solname,rhsname,bndname,modename, &
        & pcname,arg

    ! computational mesh
    type (tri_mesh) :: mesh

    ! stiffness and mass matrices
    class (sparse_matrix), allocatable :: A, B, R

    ! rhs/solution vectors
    real(kind(1d0)), allocatable :: u(:),f(:),g(:),z(:)

    ! solvers
    class(iterative_solver), allocatable :: krylov
    class(preconditioner), allocatable :: pc

    ! some other locals
    integer :: i,mc
    integer, allocatable :: p(:),q(:)

    ! variables for reading/writing netcdf
    integer :: rcode,ncid,nodesid,fid,gid,uid



!--------------------------------------------------------------------------!
! Read in mesh data                                                        !
!--------------------------------------------------------------------------!
    ! Parse the command line arguments

    rhsname = "none                            "
    bndname = "none                            "
    pcname = "none                            "

    do i=1,iargc()
        call getarg(i,arg)

        select case(trim(arg))
        case('--out')
            call getarg(i+1,solname)
        case('--mesh')
            call getarg(i+1,meshname)
        case('--rhs')
            call getarg(i+1,rhsname)
        case('--bc')
            call getarg(i+1,modename)
        case('--bnd')
            call getarg(i+1,bndname)
        case('--pc')
            call getarg(i+1,pcname)
        case('--help')
            print *, 'This program solves boundary value problems for the  '
            print *, 'Laplace operator with either Dirichlet or Robin      '
            print *, 'boundary conditions.                                 '
            print *, 'Arguments:                                           '
            print *, '   --out <path to output file>                       '
            print *, '   --mesh <path to mesh>                             '
            print *, '   --rhs <path to right-hand side>                   '
            print *, '   --bnd <path to boundary data>                     '
            print *, '   --bc <type of bc>, robin or dirichlet             '
            print *, '   --pc <type of preconditioner>, see docs for list  '
            call exit(0)
        end select
    enddo
           
    ! Read the mesh
    call read_mesh(meshname,mesh)

    allocate(csr_matrix::A)
    allocate(csr_matrix::B)
    allocate(csr_matrix::R)

    call assemble(mesh,A)
    call stiffness_matrix(mesh,A,1.d0)

    call assemble(mesh,B)
    call mass_matrix(mesh,B)



!--------------------------------------------------------------------------!
! Load in the right-hand side and boundary data                            !
!--------------------------------------------------------------------------!
    allocate( u(mesh%nn), f(mesh%nn), g(mesh%nn), z(mesh%nn) )
    u = 0.d0
    f = 0.d0
    g = 0.d0

    if (trim(rhsname) /= "none") then
        ! Get the values of the right-hand side from a file
        rcode = nf90_open(rhsname,nf90_nowrite,ncid)
        rcode = nf90_inq_varid(ncid,'u',fid)
        rcode = nf90_get_var(ncid,fid,f)
        rcode = nf90_close(ncid)

        call B%matvec(f,z)
        f = z
    endif

    if (trim(bndname) /= "none") then
        rcode = nf90_open(bndname,nf90_nowrite,ncid)
        rcode = nf90_inq_varid(ncid,'u',gid)
        rcode = nf90_get_var(ncid,gid,u)
        rcode = nf90_close(ncid)

    endif

    call A%write_to_file("a")

!--------------------------------------------------------------------------!
! Do some tests n stuff                                                    !
!--------------------------------------------------------------------------!

    allocate(p(A%nrow),q(A%nrow))

    call greedy_multicolor(A,p,mc)

!    do i=1,A%nrow
!        q( p(i) ) = i
!    enddo

    call A%permute(p)
    call A%write_to_file("ap")

end program poisson

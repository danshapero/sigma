!> @mainpage BVP program
!> @author
!> Daniel Shapero, University of Washington
!> @brief
!> Solve the Poisson problem
!> \f$ -\nabla^2u = f, \qquad u|_{\partial\Omega} = g\f$

program bvp

    use mesh_mod
    use abstract_matrix_mod
    use csr_matrix_mod
    use linalg_mod
    use fem_mod
    use netcdf

    implicit none

    ! command line arguments
    character(len=32) :: meshname,solname,rhsname,bndname,modename

    ! computational mesh
    type (tri_mesh) :: mesh

    ! stiffness and mass matrices
    type (csr_matrix) :: A, B, R

    ! rhs/solution vectors
    real(kind(1d0)), allocatable :: u(:),f(:),g(:),z(:)

    ! some other locals
    integer :: i,j,k,n

    ! variables for reading/writing netcdf
    integer :: rcode,ncid,nodesid,fid,gid,uid

!--------------------------------------------------------------------------!
! Read in mesh data                                                        !
!--------------------------------------------------------------------------!
    ! Get the command line arguments
    call getarg(1,meshname)
    call getarg(2,solname)
    call getarg(3,rhsname)
    call getarg(4,bndname)
    call getarg(5,modename)

    ! Read the mesh
    call read_mesh(meshname,mesh)

    call assemble(mesh,A)
    call stiffness_matrix(mesh,A,1.d0)

    call assemble(mesh,B)
    call mass_matrix(mesh,B)

    call A%write_to_file("a")



!--------------------------------------------------------------------------!
! Load in the right-hand side and boundary data                            !
!--------------------------------------------------------------------------!
    allocate( u(mesh%nn), f(mesh%nn), g(mesh%nn), z(mesh%nn) )
    u = 0.d0
    f = 0.d0
    g = 0.d0

    if (trim(rhsname) /= "none" ) then
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



!--------------------------------------------------------------------------!
! Solve for u using the conjugate gradient method                          !
!--------------------------------------------------------------------------!
    call allocate_linalg_workarrays( mesh%nn )

    do i=1,A%nrow
        if ( mesh%bnd(i)/=0 ) f(i) = 0.d0
    enddo

    call subset_cg(A,u,f,1.0D-8,mesh%bnd,0)



!--------------------------------------------------------------------------!
! Write u to a netcdf file                                                 !
!--------------------------------------------------------------------------!
    rcode = nf90_create(trim(solname)//'.nc',nf90_clobber,ncid)
    rcode = nf90_def_dim(ncid,'nodes',mesh%nn,nodesid)
    rcode = nf90_def_var(ncid,'u',nf90_double,nodesid,uid)
    rcode = nf90_enddef(ncid)
    rcode = nf90_close(ncid)
    rcode = nf90_open(trim(solname)//'.nc',nf90_write,ncid)
    rcode = nf90_put_var(ncid,uid,u)
    rcode = nf90_close(ncid)


end program bvp

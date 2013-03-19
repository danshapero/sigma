program poisson

    use mesh_mod
    use matrix_mod
    use fem_mod
    use linalg_mod
    use netcdf

    implicit none

    ! command line arguments
    character(len=32) :: meshname,solname,rhsname,bndname,modename

    ! computational mesh
    type (tri_mesh) :: mesh

    ! stiffness and mass matrices
    type (csr_sparse_matrix) :: A, B, R

    ! rhs/solution vectors
    real(kind=8), dimension(:), allocatable :: u,f,g,z

    ! some other locals    
    integer :: i,j,k,n,ierr

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

    print *, 'Hi there!'

    ! Assemble and fill the matrices
    call init_matrix(A, mesh%nn, mesh%nn, mesh%nn+2*mesh%nl )
    call assemble(mesh,A)
    call stiffness_matrix(mesh,1.d0,A)
    call init_matrix(B, mesh%nn, mesh%nn, mesh%nn+2*mesh%nl )
    B%ia = A%ia
    B%ja = A%ja
    call mass_matrix(mesh,B)

    print *, 'Hello again!'

    if (trim(modename) == "robin") then
        call assemble_boundary(mesh,R)
        call robin_matrix(mesh,R)
        call sub_matrix_add(A,R,A,ierr)
    endif

    print *, 'Salutations!'



!--------------------------------------------------------------------------!
! Load in the right-hand side and boundary data                            !
!--------------------------------------------------------------------------!
    allocate( u(mesh%nn), f(mesh%nn), g(mesh%nn), z(mesh%nn) )
    u = 0.d0
    f = 0.d0
    g = 0.d0

    print *, 'Bonjour!'

    if (trim(rhsname) /= "none") then
        ! Get the values of the right-hand side from a file
        rcode = nf90_open(rhsname,nf90_nowrite,ncid)
        rcode = nf90_inq_varid(ncid,'u',fid)
        rcode = nf90_get_var(ncid,fid,f)
        rcode = nf90_close(ncid)

        call B%matvec(f,z)
        f = z
    endif

    print *, 'Sup!'

    if (trim(bndname) /= "none") then
        rcode = nf90_open(bndname,nf90_nowrite,ncid)
        rcode = nf90_inq_varid(ncid,'u',gid)
        rcode = nf90_get_var(ncid,gid,u)
        rcode = nf90_close(ncid)

        if (trim(modename) == "robin") then
            call R%matvec(u,z)
            f = f+z
        endif
    endif

    print *, 'Yo!'

!--------------------------------------------------------------------------!
! Solve for u using the conjugate gradient method                          !
!--------------------------------------------------------------------------!
    call allocate_linalg_workarrays( mesh%nn )

    if (trim(modename) == "dirichlet") then

        print *, minval(f),maxval(f)

        do n=1,A%nrow
            if ( mesh%bnd(n) /= 0 ) f(n) = 0.d0
        enddo

        print *, 'In Xanadu'

        call subset_cg(A,u,f,1.0E-8,mesh%bnd,0)
    elseif (trim(modename) == "robin") then
        call cg(A,u,f,1.0E-8)
    endif

    print *, minval(u),maxval(u)

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







end program poisson

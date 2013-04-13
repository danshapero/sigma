!--------------------------------------------------------------------------!
module fem_mod                                                             !
!--------------------------------------------------------------------------!
    use mesh_mod
    use sparse_matrix_mod

    implicit none

contains


!--------------------------------------------------------------------------!
subroutine assemble(mesh,A)                                                !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type (tri_mesh), intent(in) :: mesh
    class (sparse_matrix), intent(inout) :: A
    ! local variables
    integer :: i,lastindex,nn,nl
    integer, dimension(2) :: edge
    integer, dimension( mesh%nn+2*mesh%nl ) :: rows,cols

    nn = mesh%nn
    nl = mesh%nl

    do i=1,nn
        rows(i) = i
        cols(i) = i
    enddo

    lastindex = nn+1

    do i=1,nl
        edge = mesh%edge(:,i)
        rows(lastindex) = edge(1)
        cols(lastindex) = edge(2)

        rows(lastindex+1) = edge(2)
        cols(lastindex+1) = edge(1)

        lastindex = lastindex+2
    enddo

    call A%init_matrix(nn,nn,nn+2*nl,rows,cols)

end subroutine assemble



!--------------------------------------------------------------------------!
subroutine stiffness_matrix(mesh,A,kappa)                                  !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type (tri_mesh), intent(in) :: mesh
    real(kind(1d0)), intent(in) :: kappa
    class (sparse_matrix), intent(inout) :: A
    ! local variables
    integer :: n
    integer, dimension(3) :: elem
    real(kind(1d0)) :: area
    real(kind(1d0)), dimension(3,3) :: AE
    real(kind(1d0)), dimension(3,2) :: D

    do n=1,mesh%ne
        !--------------------------------------
        ! Compute the element stiffness matrix
        elem = mesh%elem(:,n)
        D(1,:) = mesh%x(:,elem(3))-mesh%x(:,elem(2))
        D(2,:) = mesh%x(:,elem(1))-mesh%x(:,elem(3))
        D(3,:) = mesh%x(:,elem(2))-mesh%x(:,elem(1))
        area = 0.5*dabs( D(1,1)*D(2,2)-D(1,2)*D(2,1) )
        AE = 0.25/area*kappa*matmul(D,transpose(D))
        ! See "Numerical Treatment of Partial Differential Equations" by
        ! Grossman, Roos and Stynes for the derivation of this formula for
        ! the element stiffness matrix.

        call A%add_values(elem,elem,AE)

    enddo

    A%pos_def = .true.
    A%symmetric = .true.
    A%diag_dominant = .true.

end subroutine stiffness_matrix



!--------------------------------------------------------------------------!
subroutine mass_matrix(mesh,B)                                             !
!--------------------------------------------------------------------------!
! Fill in the entries of the sparse matrix B constituting the mass matrix  !
!   for C0 finite elements on the tri_mesh mesh                            !
! Input: a tri_mesh mesh                                                   !
! Output: the sparse_matrix B, with entries filled                         !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type (tri_mesh), intent(in) :: mesh
    class (sparse_matrix), intent(inout) :: B
    ! local variables
    integer :: i,n
    integer, dimension(3) :: elem
    real(kind(1d0)) :: area
    real(kind(1d0)), dimension(3,3) :: BE

    do n=1,mesh%ne
        !---------------------------------
        ! Compute the element mass matrix
        elem = mesh%elem(:,n)
        BE(1,1:2) = mesh%x(:,elem(2))-mesh%x(:,elem(3))
        BE(2,1:2) = mesh%x(:,elem(1))-mesh%x(:,elem(3))
        area = 0.5*dabs( BE(1,1)*BE(2,2)-BE(1,2)*BE(2,1) )
        BE = area/12
        do i=1,3
            BE(i,i) = BE(i,i)+area/12
        enddo

        call B%add_values(elem,elem,BE)

    enddo

    B%pos_def = .true.
    B%symmetric = .true.
    B%diag_dominant = .true.

end subroutine mass_matrix



!--------------------------------------------------------------------------!
function matvec_mass_matrix(mesh,u)                                        !
!--------------------------------------------------------------------------!
! Multiply a vector u by the mass matrix for a tri_mesh mesh.              !
! This is useful when we only need one multiplication by the mass matrix,  !
!   e.g. to compute the load vector for the FEM, rather than computing and !
!   storing the whole mass matrix.                                         !
! Input: a tri_mesh mesh                                                   !
!        a vector u                                                        !
! Output: a vector y, the product of the stiffness matrix and u            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type (tri_mesh), intent(in) :: mesh
    real(kind(1d0)), dimension( mesh%nn ), intent(in) :: u
    real(kind(1d0)), dimension( mesh%nn ) :: matvec_mass_matrix
    ! local variables
    integer :: i,n
    integer, dimension(3) :: elem
    real(kind(1d0)) :: area
    real(kind(1d0)), dimension(3) :: y
    real(kind(1d0)), dimension(3,3) :: BE

    matvec_mass_matrix = 0.d0
    do n=1,mesh%ne
        !---------------------------------
        ! Compute the element mass matrix
        elem = mesh%elem(:,n)
        BE(1,1:2) = mesh%x(:,elem(2))-mesh%x(:,elem(3))
        BE(2,1:2) = mesh%x(:,elem(1))-mesh%x(:,elem(3))
        area = 0.5*dabs( BE(1,1)*BE(2,2)-BE(1,2)*BE(2,1) )
        BE = area/12
        do i=1,3
            BE(i,i) = BE(i,i)+area/12
        enddo

        !--------------------------------------------------
        ! Multiply the components of u from the current
        ! element by the element mass matrix and them to y
        y = matmul( BE,u(elem) )
        matvec_mass_matrix(elem) = matvec_mass_matrix(elem)+y
    enddo

end function matvec_mass_matrix



!--------------------------------------------------------------------------!
subroutine assemble_boundary(mesh,R)                                       !
!--------------------------------------------------------------------------!
! Fill in the arrays ia,ja of a sparse_matrix R which describe its non-    !
!   zero structure for the tri_mesh mesh, provided R only operatoes on     !
!   boundary nodes; for example, R is the matrix discretizing Robin        !
!   boundary conditions.                                                   !
! Input: a tri_mesh mesh                                                   !
!        an allocated sparse_matrix R                                      !
! Output: the sparse_matrix R, with its arrays ia,ja properly filled       !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type (tri_mesh), intent(in) :: mesh
    class (sparse_matrix), intent(inout) :: R
    ! local variables
    integer :: i,lastindex,nnzr
    integer, dimension(2) :: edge
    integer, allocatable :: rows(:),cols(:)


    !-------------------------------------------------------------
    ! #non-zero entries in R = #boundary nodes + 2*boundary edges
    nnzr = sum(abs(mesh%bnd))+2*sum(abs(mesh%bnd_edge))
    allocate(rows(nnzr),cols(nnzr))

    lastindex = 1
    do i=1,mesh%nn
        if ( mesh%bnd(i)/=0 ) then
            rows(lastindex) = i
            cols(lastindex) = i
            lastindex = lastindex+1
        endif
    enddo

    do i=1,mesh%nl
        edge = mesh%edge(:,i)
        if (mesh%bnd_edge(i) /= 0) then
            rows( lastindex ) = edge(1)
            cols( lastindex ) = edge(2)

            rows( lastindex+1 ) = edge(2)
            cols( lastindex+1 ) = edge(1)

            lastindex = lastindex+2
        endif
    enddo

    call R%init_matrix( mesh%nn, mesh%nn, nnzr, rows, cols )
    deallocate(rows,cols)
    

end subroutine assemble_boundary



!--------------------------------------------------------------------------!
subroutine robin_matrix(mesh,R)                                            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type (tri_mesh), intent(in) :: mesh
    class (sparse_matrix), intent(inout) :: R
    ! local variables
    integer :: i
    integer, dimension(2) :: edge
    real(kind(1d0)) :: dx
    real(kind(1d0)), dimension(2,2) :: z,RE

    do i=1,mesh%nl
        if (mesh%bnd_edge(i) /= 0) then
            edge = mesh%edge(:,i)
            z(:,1) = mesh%x(:,edge(1))
            z(:,2) = mesh%x(:,edge(2))
            dx = dsqrt( (z(1,1)-z(1,2))**2 + (z(2,1)-z(2,2))**2 )
            RE = dx/6
            RE(1,1) = dx/3
            RE(2,2) = dx/3

            call R%add_values(edge,edge,RE)
        endif
    enddo

    R%pos_def = .true.
    R%symmetric = .true.
    R%diag_dominant = .true.

end subroutine robin_matrix



!--------------------------------------------------------------------------!
function gradient(mesh,u,n)                                                !
!--------------------------------------------------------------------------!
! Compute the gradient  of a function u defined on the tri_mesh mesh.      !
! Input: a tri_mesh mesh                                                   !
!        a vector u                                                        !
!        an integer n                                                      !
! Output: a vector gradient(2), giving the gradient of u within element n  !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type (tri_mesh), intent(in) :: mesh
    real(kind(1d0)), dimension(mesh%nn), intent(in) :: u
    integer, intent(in) :: n
    real(kind(1d0)), dimension(2) :: gradient
    ! local variables
    integer, dimension(3) :: elem
    real(kind(1d0)) :: dx,det
    real(kind(1d0)), dimension(2,2) :: S,Q

    elem = mesh%elem(:,n)

    ! z -> S*z+x_3 maps the reference triangle to the physical triangle
    S(:,1) = mesh%x(:,elem(1))-mesh%x(:,elem(3))
    S(:,2) = mesh%x(:,elem(2))-mesh%x(:,elem(3))

    ! x -> Q*(x-x_3) maps the physical triangle to the reference triangle
    det = S(1,1)*S(2,2)-S(1,2)*S(2,1)
    Q(1,1) = S(2,2)/det
    Q(2,2) = S(1,1)/det
    Q(2,1) = -S(2,1)/det
    Q(1,2) = -S(1,2)/det

    ! Compute the directional derivatives of u along the edges of the
    ! physical triangle
    dx = dsqrt( S(1,1)**2 + S(2,1)**2 )
    gradient(1) = ( u(elem(1))-u(elem(3)) )/dx
    dx = dsqrt( S(1,2)**2 + S(2,2)**2 )
    gradient(2) = ( u(elem(2))-u(elem(3)) )/dx

    ! Multiply by Q to find the gradient vector
    gradient = matmul(gradient,Q)
    

end function gradient



end module fem_mod

module matrix_c

    use iso_c_binding
    use ellpack
    use bsr
    use csr

    implicit none


!--------------------------------------------------------------------------!
! C wrapper to Fortran matrix data types                                   !
!--------------------------------------------------------------------------!
type, bind(c) :: sparse_matrix_c
    type(c_ptr) :: p
    integer(c_int) :: nrow,ncol,mat_type
end type sparse_matrix_c



contains






!==========================================================================!
!==========================================================================!
!==== Subroutines to get Fortran pointers from C pointers to matrices  ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine get_sparse_matrix_c(cmat,mat_type) bind(c)                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(sparse_matrix_c), intent(inout) :: cmat
    integer(c_int), intent(in), value :: mat_type
    ! local variables
    type(csr_matrix), pointer :: csr_mat
    type(bsr_matrix), pointer :: bsr_mat
    type(ellpack_matrix), pointer :: ellpack_mat

    cmat%mat_type = mat_type
    select case(mat_type)
        case(0)
            allocate(csr_mat)
            cmat%p = c_loc(csr_mat)
        case(1)
            allocate(bsr_mat)
            cmat%p = c_loc(bsr_mat)
        case(2)
            allocate(ellpack_mat)
            cmat%p = c_loc(ellpack_mat)
    end select

end subroutine get_sparse_matrix_c



!--------------------------------------------------------------------------!
subroutine sparse_matrix_f(fmat,cmat)                                      !
!--------------------------------------------------------------------------!
    implicit none
    class(sparse_matrix), pointer, intent(inout) :: fmat
    type(sparse_matrix_c), intent(in) :: cmat
    ! local variables
    type(csr_matrix), pointer :: csr_mat
    type(bsr_matrix), pointer :: bsr_mat
    type(ellpack_matrix), pointer :: ellpack_mat

    select case(cmat%mat_type)
        case(0)
            allocate(csr_mat)
            call c_f_pointer(cmat%p,csr_mat)
            fmat => csr_mat
        case(1)
            allocate(bsr_mat)
            call c_f_pointer(cmat%p,bsr_mat)
            fmat => bsr_mat
        case(2)
            allocate(ellpack_mat)
            call c_f_pointer(cmat%p,ellpack_mat)
            fmat => ellpack_mat
    end select

end subroutine sparse_matrix_f






!==========================================================================!
!==========================================================================!
!==== C wrapper subroutines to sparse matrix type-bound procedures     ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine init_c(cmat,nrow,ncol,nnz) bind(c)                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(sparse_matrix_c), intent(inout) :: cmat
    integer(c_int), intent(in), value :: nrow,ncol,nnz
    ! local variables
    class(sparse_matrix), pointer :: A

    call sparse_matrix_f(A,cmat)
    call A%init(nrow,ncol,nnz)

end subroutine init_c



!--------------------------------------------------------------------------!
subroutine build_c(cmat,rows,cols,nnz) bind(c)                             !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(sparse_matrix_c), intent(in) :: cmat
    integer(c_int), intent(in) :: rows(nnz), cols(nnz)
    integer(c_int), intent(in), value :: nnz
    ! local variables
    class(sparse_matrix), pointer :: A

    call sparse_matrix_f(A,cmat)
    call A%build(rows+1,cols+1)

end subroutine build_c



!--------------------------------------------------------------------------!
subroutine get_value_c(cmat,i,j,val) bind(c)                               !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(sparse_matrix_c), intent(in) :: cmat
    integer(c_int), intent(in), value :: i,j
    real(c_double), intent(out) :: val
    ! local variables
    class(sparse_matrix), pointer :: A

    call sparse_matrix_f(A,cmat)
    val = A%get_value(i+1,j+1)

end subroutine get_value_c



!--------------------------------------------------------------------------!
subroutine set_value_c(cmat,i,j,val) bind(c)                               !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(sparse_matrix_c), intent(in) :: cmat
    integer(c_int), intent(in), value :: i,j
    real(c_double), intent(in), value :: val
    ! local variables
    class(sparse_matrix), pointer :: A

    call sparse_matrix_f(A,cmat)
    call A%set_value(i+1,j+1,val)

end subroutine set_value_c




end module matrix_c

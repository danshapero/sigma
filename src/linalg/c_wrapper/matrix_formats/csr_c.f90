module csr_c

    use iso_c_binding
    use csr

    implicit none

contains



!--------------------------------------------------------------------------!
subroutine csr_matrix_c(ptr) bind(c)                                       !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(c_ptr), intent(out) :: ptr
    ! local variables
    type(csr_matrix), pointer :: A

    allocate(A)
    ptr = c_loc(A)

end subroutine csr_matrix_c



!--------------------------------------------------------------------------!
subroutine csr_init_c(ptr,nrow,ncol,nnz) bind(c)                           !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in), value :: nrow,ncol,nnz
    ! local variables
    type(csr_matrix), pointer :: A

    allocate(A)
    call c_f_pointer(ptr,A)
    call A%init(nrow,ncol,nnz)
    
end subroutine csr_init_c



!--------------------------------------------------------------------------!
subroutine csr_build_c(ptr,rows,cols,nnz) bind(c)                          !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(c_ptr), intent(inout) :: ptr
    integer(c_int), intent(in) :: rows(nnz),cols(nnz)
    integer(c_int), intent(in), value :: nnz
    ! local variables
    type(csr_matrix), pointer :: A

    allocate(A)
    call c_f_pointer(ptr,A)
    call A%build(rows+1,cols+1)

end subroutine csr_build_c



!--------------------------------------------------------------------------!
subroutine csr_get_value_c(ptr,i,j,val) bind(c)                            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in), value :: i,j
    real(c_double), intent(out) :: val
    ! local variables
    type(csr_matrix), pointer :: A

    allocate(A)
    call c_f_pointer(ptr,A)
    val = A%get_value(i+1,j+1)

end subroutine csr_get_value_c



!--------------------------------------------------------------------------!
subroutine csr_set_value_c(ptr,i,j,val) bind(c)                            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in), value :: i,j
    real(c_double), intent(in), value :: val
    ! local variables
    type(csr_matrix), pointer :: A

    allocate(A)
    call c_f_pointer(ptr,A)
    call A%set_value(i+1,j+1,val)

end subroutine csr_set_value_c



!--------------------------------------------------------------------------!
subroutine csr_matvec_c(ptr,x,y,m,n) bind(c)                               !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(out) :: y(m)
    integer(c_int), intent(in), value :: m,n
    ! local variables
    type(csr_matrix), pointer :: A

    allocate(A)
    call c_f_pointer(ptr,A)
    call A%matvec(x,y)

end subroutine csr_matvec_c



end module csr_c

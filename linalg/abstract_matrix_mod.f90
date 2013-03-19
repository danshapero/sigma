module abstract_matrix_mod

	implicit none



type, abstract :: sparse_matrix
    integer :: nrow,ncol,nnz,max_degree
    logical :: symmetric,pos_def,m_matrix,diagonally_dominant
contains
    procedure(init_matrix_interface), deferred :: init_matrix
    procedure(get_value_interface), deferred :: get_value
    procedure(get_values_interface), deferred :: get_values
    procedure(set_value_interface), deferred :: set_value, add_value
    procedure(set_values_interface), deferred :: set_values, add_values
    procedure(matvec_interface), deferred :: matvec
    procedure(subblock_matvec_interface), deferred :: subblock_matvec
    procedure(subset_matvec_interface), deferred :: subset_matvec
    procedure(write_to_file_interface), deferred :: write_to_file
end type sparse_matrix



abstract interface

subroutine init_matrix_interface(A,nrow,ncol,nnz,rows,cols)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: nrow,ncol,nnz
    integer, intent(in), optional :: rows(:), cols(:)
end subroutine init_matrix_interface


real(kind(1d0)) function get_value_interface(A,i,j) result(value)
    import :: sparse_matrix
    class (sparse_matrix), intent(in) :: A
    integer, intent(in) :: i,j
end function get_value_interface


function get_values_interface(A,rows,cols)
    import :: sparse_matrix
    class (sparse_matrix), intent(in) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)) :: get_values_interface(size(rows),size(cols))
end function get_values_interface


subroutine set_value_interface(A,i,j,value)
    import :: sparse_matrix
    class (sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: value
end subroutine set_value_interface


subroutine set_values_interface(A,rows,cols,values)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: values(size(rows),size(cols))
end subroutine set_values_interface


subroutine matvec_interface(A,x,y)
    import :: sparse_matrix
    class (sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(in) :: x(:)
    real(kind(1d0)), intent(out) :: y(:)
end subroutine matvec_interface


subroutine subblock_matvec_interface(A,x,y,i1,i2,j1,j2)
    import :: sparse_matrix
    class (sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(in) :: x(:)
    real(kind(1d0)), intent(inout) :: y(:)
    integer, intent(in) :: i1,i2,j1,j2
end subroutine subblock_matvec_interface


subroutine subset_matvec_interface(A,x,y,setlist,set1,set2)
    import :: sparse_matrix
    class (sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(in) :: x(:)
    real(kind(1d0)), intent(inout) :: y(:)
    integer, intent(in) :: setlist(:), set1, set2
end subroutine subset_matvec_interface


subroutine write_to_file_interface(A,filename)
    import :: sparse_matrix
    class (sparse_matrix), intent(in) :: A
    character(len=*), intent(in) :: filename
end subroutine write_to_file_interface


end interface




contains

!-----------
! Nothing! This module is just for defining the abstract sparse_matrix type


end module abstract_matrix_mod

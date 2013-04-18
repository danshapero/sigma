module sparse_matrix_mod

    implicit none



type, abstract :: sparse_matrix
    integer :: nrow,ncol,nnz,max_degree
    logical :: symmetric,pos_def,m_matrix,diag_dominant
contains
    ! Constructor and accessors/mutators
    procedure(init_matrix_interface), deferred :: init_matrix
    procedure(get_value_interface), deferred :: get_value
    procedure(get_values_interface), deferred :: get_values
    procedure(get_neighbors_interface), deferred :: get_neighbors
    procedure(set_value_interface), deferred :: set_value, add_value
    procedure(set_values_interface), deferred :: set_values, add_values
    procedure(permute_interface), deferred :: permute
    procedure(subset_matrix_add_interface), deferred :: subset_matrix_add
    ! matrix multiplication routines
    procedure(matvec_interface), deferred :: matvec
    ! forward- and back-solves for triangular systems
    procedure(trisolve_interface), deferred :: backsolve, forwardsolve
    ! routines for i/o and validation
    procedure(convert_to_coo_interface), deferred :: convert_to_coo
    procedure(write_to_file_interface), deferred :: write_to_file
end type sparse_matrix



abstract interface

subroutine init_matrix_interface(A,nrow,ncol,nnz,rows,cols,vals)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: nrow,ncol,nnz
    integer, intent(in), optional :: rows(:),cols(:)
    real(kind(1d0)), intent(in), optional :: vals(:)
end subroutine init_matrix_interface


real(kind(1d0)) function get_value_interface(A,i,j) result(value)
    import :: sparse_matrix
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: i,j
end function get_value_interface


function get_values_interface(A,rows,cols)
    import :: sparse_matrix
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)) :: get_values_interface(size(rows),size(cols))
end function get_values_interface


function get_neighbors_interface(A,row)
    import :: sparse_matrix
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: row
    integer :: get_neighbors_interface( A%max_degree )
end function get_neighbors_interface


subroutine set_value_interface(A,i,j,value)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: value
end subroutine set_value_interface


subroutine set_values_interface(A,rows,cols,values)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: values(size(rows),size(cols))
end subroutine set_values_interface


subroutine permute_interface(A,p)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
end subroutine permute_interface


subroutine subset_matrix_add_interface(A,B)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    class(sparse_matrix), intent(in) :: B
end subroutine


subroutine matvec_interface(A,x,y,rows,cols)
    import :: sparse_matrix
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(in) :: x(:)
    real(kind(1d0)), intent(out) :: y(:)
    integer, intent(in), optional :: rows(2),cols(2)
end subroutine matvec_interface


subroutine trisolve_interface(A,x)
    import :: sparse_matrix
    class(sparse_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
end subroutine


subroutine convert_to_coo_interface(A,rows,cols,vals)
    import :: sparse_matrix
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: rows(:),cols(:)
    real(kind(1d0)), intent(out), optional :: vals(:)
end subroutine


subroutine write_to_file_interface(A,filename)
    import :: sparse_matrix
    class(sparse_matrix), intent(in) :: A
    character(len=*), intent(in) :: filename
end subroutine write_to_file_interface


end interface




contains

!-----------
! Nothing! This module is just for defining the abstract sparse_matrix type


end module sparse_matrix_mod

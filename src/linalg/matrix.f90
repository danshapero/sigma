module matrix

    implicit none



type, abstract :: sparse_matrix
    integer :: nrow,ncol,nnz,max_degree
    logical :: symmetric,pos_def,m_matrix,diag_dominant
contains
    ! Constructor and accessors/mutators
	procedure(init_interface), deferred :: init
    procedure(build_interface), deferred :: build
    procedure(get_value_interface), deferred :: get_value
    procedure(get_values_interface), deferred :: get_values
    procedure(get_neighbors_interface), deferred :: get_neighbors
    procedure(set_value_interface), deferred :: set_value, add_value
    procedure(set_values_interface), deferred :: set_values, add_values
    procedure(zero_interface), deferred :: zero
    procedure(permute_interface), deferred :: permute
    procedure(subset_matrix_add_interface), deferred :: subset_matrix_add
    ! matrix multiplication routines
    procedure(matvec_interface), deferred :: matvec
    ! forward- and back-solves for triangular systems
    procedure(trisolve_interface), deferred :: backsolve, forwardsolve
    ! routines for i/o and validation
    procedure(convert_to_coo_interface), deferred :: convert_to_coo
    procedure :: write_to_file
end type sparse_matrix



abstract interface

subroutine init_interface(A,nrow,ncol,nnz,rows,cols,params)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: nrow,ncol,nnz
    integer, intent(in), optional :: rows(:),cols(:),params(:)
end subroutine init_interface


subroutine build_interface(A,rows,cols)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
end subroutine build_interface


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


subroutine set_value_interface(A,i,j,val)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: val
end subroutine set_value_interface


subroutine set_values_interface(A,rows,cols,vals)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: vals(size(rows),size(cols))
end subroutine set_values_interface


subroutine zero_interface(A)
    import :: sparse_matrix
    class(sparse_matrix), intent(inout) :: A
end subroutine zero_interface


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



end interface




contains



!--------------------------------------------------------------------------!
subroutine write_to_file(A,filename)                                       !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    character(len=*), intent(in) :: filename
    ! local variables
    integer :: i,rows(A%nnz),cols(A%nnz)
    real(kind(1d0)) :: vals(A%nnz)

    call A%convert_to_coo(rows,cols,vals)
    open(unit=100,file=trim(filename)//'.txt')
    write(100,*) A%nrow,A%ncol,A%nnz
    do i=1,A%nnz
        write(100,*) rows(i),cols(i),vals(i)
    enddo
    close(100)

end subroutine write_to_file





end module matrix

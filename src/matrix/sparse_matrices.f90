!==========================================================================!
!==========================================================================!
module sparse_matrices                                                     !
!==========================================================================!
!==========================================================================!

use graphs
use sparse_matrix_interfaces
use default_matrices
use cs_matrices
use ellpack_matrices
use sparse_matrix_algebra, only: sparse_matrix_sum, sparse_matrix_product

implicit none


! Parameter giving the number of different matrix formats available.
! Test programs iterate over every matrix type to check correctness of each
! matrix implementation.
integer, parameter :: num_matrix_types = 5



interface choose_matrix_type
    module procedure choose_matrix_type_by_int, choose_matrix_type_by_name
end interface

contains



!==========================================================================!
!==== Matrix factory methods                                           ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine choose_matrix_type_by_int(A, t)                                 !
!--------------------------------------------------------------------------!
    class(sparse_matrix_interface), pointer, intent(inout) :: A
    integer, intent(in) :: t

    nullify(A)

    select case(t)
        case(1)
            allocate(default_row_matrix :: A)

        case(2)
            allocate(default_column_matrix :: A)

        case(3)
            allocate(csr_matrix :: A)

        case(4)
            allocate(csc_matrix :: A)

        case(5)
            allocate(ellpack_matrix :: A)

    end select

end subroutine choose_matrix_type_by_int



!--------------------------------------------------------------------------!
subroutine choose_matrix_type_by_name(A, frmt)                             !
!--------------------------------------------------------------------------!
    class(sparse_matrix_interface), pointer, intent(inout) :: A
    character(len=*), intent(in) :: frmt

    nullify(A)

    select case( trim(frmt) )
        case("default row")
            allocate(default_row_matrix :: A)

        case("default column", "default col")
            allocate(default_column_matrix :: A)

        case("csr", "crs", "harwell-boeing", "compressed sparse row")
            allocate(csr_matrix :: A)

        case("csc", "ccs", "compressed sparse column")
            allocate(csc_matrix :: A)

        case("ellpack", "itpack")
            allocate(ellpack_matrix :: A)
    end select

end subroutine choose_matrix_type_by_name



end module sparse_matrices


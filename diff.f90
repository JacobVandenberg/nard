! Created by  on 6/7/21.

module diff
    use sparse_matrices
    implicit none

    contains
        subroutine tridiag(a, b, c, T)
            ! retuns T = NxN tridiagonal matrix with a on main diagonal,
            ! b on upper diagonal and c on lower diagonal.
            !
            ! inputs: a: value on main diagonal
            !         b: value on upper diagonal
            !         c: value on lower diagonal
            ! output:
            !         T: returned matrix
            implicit none
            real (kind=8), intent(in) :: a, b, c
            real (kind=8), dimension(:, :), intent(out) :: T
            integer :: N, i
            N = size(T, 1)

            if (N==0) then
                return
            end if
            if (N==1) then
                T(1, 1) = a
                return
            end if

            ! deal with edges
            T(1, 1) = a
            T(N, N) = a
            T(1, 2) = b
            T(N, N-1) = c

            ! deal with the rest
            do i = 2, N-1
                T(i, i) = a
                T(i, i+1) = b
                T(i, i-1) = c
            end do
            return
        end subroutine tridiag

        subroutine sparse_tridiag(a, b, c, N, val, indx, jndx, ierr)
            ! retuns NxN tridiagonal sparse matrix with a on main diagonal,
            ! b on upper diagonal and c on lower diagonal.
            !
            ! this function allocates the memory for the sparse matrix
            ! which must be deallocated later
            !
            ! inputs: a: value on main diagonal
            !         b: value on upper diagonal
            !         c: value on lower diagonal
            !         N: size of output matrix
            !
            ! output: val: non-zero values in the matrix (length 3N - 2)
            !         indx: i indices (length 3N - 2)
            !         jndx: j indices (length 3N - 2)
            !         ierr: error flag (length 3N - 2)
            !

            implicit none
            real (kind=8), intent(in) :: a, b, c
            integer, intent(in) :: N
            integer, intent(out) :: ierr
            real (kind=8), dimension(:), allocatable :: val
            integer, dimension(:), allocatable :: indx, jndx
            integer :: i

            allocate(indx(3*N-2), jndx(3*N-2), val(3*N-2), STAT=ierr)
            if (ierr/=0) then
                return
            end if

            if (N == 0) then
                return
            end if

            if (N == 1) then
                val(1) = a
                indx(1) = 1
                jndx(1) = 1
                return
            end if

            val(1) = a
            val(2) = b
            indx(1) = 1
            jndx(1) = 1
            indx(2) = 1
            jndx(2) = 2
            do i = 2, N-1
                indx(3*i-3:3*i-1) = i
                jndx(3*i-3:3*i-1) = (/ i-1, i, i+1 /)
                val(3*i-3:3*i-1) = (/ c, a, b /)
            end do
            val(3*N-2) = a
            val(3*N-3) = c
            indx(3*N-2) = N
            jndx(3*N-2) = N
            indx(3*N-3) = N
            jndx(3*N-3) = N-1
            return
        end subroutine sparse_tridiag

        type (coo_matrix) function double_diff_FD(x, BC, ierr)
            ! returns a sparse double derivative matrix
            ! input:
            !   x: equispaced grid
            !   BC: integer, 1 if boundary conditions are periodic
            ! return: coo_matrix of the double derivative matrix.
            !
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            real (kind = 8), dimension(:), intent(in) :: x
            integer, intent(in) :: BC

            ! outputs
            integer, intent(out) :: ierr
            type (coo_matrix) :: return_matrix

            ! runtime
            real (kind = 8) :: dx
            integer :: xlen
            ! END DECLARATIONS

            ! BEGIN FUNCTION
            ! make sure youve got enough points
            xlen = size(x, 1)
            if (xlen < 2) then
                ierr = -1
                return
            end if

            dx = x(2) - x(1)
            call sparse_tridiag(dble(-2)/ (dx*dx), dble(1)/(dx*dx), dble(1)/(dx*dx), xlen,&
                    return_matrix%vals, return_matrix%indx, return_matrix%jndx, ierr)
            return_matrix%n = xlen; return_matrix%m = xlen

            if (BC == 1) then
                call return_matrix%set_value(1, xlen, return_matrix%get_value(1, 2), ierr)
                call return_matrix%set_value(xlen, 1, return_matrix%get_value(xlen, xlen-1), ierr)
            else
                call return_matrix%set_value(1, 2, -return_matrix%get_value(1, 1), ierr)
                call return_matrix%set_value(xlen, xlen-1, -return_matrix%get_value(xlen, xlen), ierr)
            end if

            double_diff_FD = return_matrix
            return
            ! END FUNCTION
        end function double_diff_FD

        type (coo_matrix) function laplacian2D(x, y, BCx, BCy, ierr)
            !
            ! generates a sparse laplacian matrix
            !
            ! inputs:
            !   x, y: vector of equispaced grid points
            !   BCx, BCy: boundary conditions for x and y directions (1 if periodic)
            !
            ! outputs:
            !   ierr: error flag
            !   return value: laplacian matrix in sparse coo form
            !
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            real (kind=8), dimension(:), intent(in) :: x, y
            integer, intent(in) :: BCx, BCy

            ! outputs
            integer, intent(out) :: ierr

            ! runtime
            type (coo_matrix) :: DDx_matrix, DDy_matrix


            ! END DECLARATIONS

            ! BEGIN FUNCTION
            DDx_matrix = double_diff_FD(x, BCx, ierr)
            DDy_matrix = double_diff_FD(y, BCy, ierr)

            DDx_matrix = sparse_kron(DDx_matrix, speye(size(y, 1), ierr), ierr)
            DDy_matrix = sparse_kron(speye(size(x, 1), ierr), DDy_matrix, ierr)
            call sparse_add(DDx_matrix, DDy_matrix, ierr)
            laplacian2D = DDx_matrix
        end function laplacian2D

end module diff
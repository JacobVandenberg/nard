! Created by  on 6/7/21.

module diff
    use sparse_matrices
    use precision
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
            real (rp), intent(in) :: a, b, c
            real (rp), dimension(:, :), intent(out) :: T
            integer (ip) :: N, i
            N = size(T, 1_ip)

            if (N==0_ip) then
                return
            end if
            if (N==1_ip) then
                T(1_ip, 1_ip) = a
                return
            end if

            ! deal with edges
            T(1_ip, 1_ip) = a
            T(N, N) = a
            T(1_ip, 2_ip) = b
            T(N, N-1_ip) = c

            ! deal with the rest
            do i = 2_ip, N-1_ip
                T(i, i) = a
                T(i, i+1_ip) = b
                T(i, i-1_ip) = c
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
            real (rp), intent(in) :: a, b, c
            integer (ip), intent(in) :: N
            integer (ip), intent(out) :: ierr
            real (rp), dimension(:), allocatable :: val
            integer (ip), dimension(:), allocatable :: indx, jndx
            integer (ip) :: i

            allocate(indx(3_ip*N-2_ip), jndx(3_ip*N-2_ip), val(3_ip*N-2_ip), STAT=ierr)
            if (ierr/=0_ip) then
                return
            end if

            if (N == 0_ip) then
                return
            end if

            if (N == 1_ip) then
                val(1_ip) = a
                indx(1_ip) = 1_ip
                jndx(1_ip) = 1_ip
                return
            end if

            val(1_ip) = a
            val(2_ip) = b
            indx(1_ip) = 1_ip
            jndx(1_ip) = 1_ip
            indx(2_ip) = 1_ip
            jndx(2_ip) = 2_ip
            do i = 2_ip, N-1_ip
                indx(3_ip*i-3_ip:3_ip*i-1_ip) = i
                jndx(3_ip*i-3_ip:3_ip*i-1_ip) = (/ i-1_ip, i, i+1_ip /)
                val(3_ip*i-3_ip:3_ip*i-1_ip) = (/ c, a, b /)
            end do
            val(3_ip*N-2_ip) = a
            val(3_ip*N-3_ip) = c
            indx(3_ip*N-2_ip) = N
            jndx(3_ip*N-2_ip) = N
            indx(3_ip*N-3_ip) = N
            jndx(3_ip*N-3_ip) = N-1_ip
            return
        end subroutine sparse_tridiag

        type (coo_matrix) function double_diff_FD(x, BC, ierr)
            ! returns a sparse double derivative matrix
            ! input:
            !   x: equispaced grid
            !   BC: integer (ip, 1 if boundary conditions are periodic
            ! return: coo_matrix of the double derivative matrix.
            !
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            real (rp), dimension(:), intent(in) :: x
            integer (ip), intent(in) :: BC

            ! outputs
            integer (ip), intent(out) :: ierr
            type (coo_matrix) :: return_matrix

            ! runtime
            real (rp) :: dx
            integer (ip) :: xlen
            ! END DECLARATIONS

            ! BEGIN FUNCTION
            ! make sure youve got enough points
            xlen = size(x, 1_ip)
            if (xlen < 2_ip) then
                ierr = -1_ip
                return
            end if

            dx = x(2_ip) - x(1_ip)
            call sparse_tridiag(dble(-2)/ (dx*dx), dble(1)/(dx*dx), dble(1)/(dx*dx), xlen,&
                    return_matrix%vals, return_matrix%indx, return_matrix%jndx, ierr)
            return_matrix%n = xlen; return_matrix%m = xlen

            if (BC == 1_ip) then
                call return_matrix%set_value(1_ip, xlen, return_matrix%get_value(1_ip, 2_ip), ierr)
                call return_matrix%set_value(xlen, 1_ip, return_matrix%get_value(xlen, xlen-1_ip), ierr)
            else
                call return_matrix%set_value(1_ip, 2_ip, -return_matrix%get_value(1_ip, 1_ip), ierr)
                call return_matrix%set_value(xlen, xlen-1_ip, -return_matrix%get_value(xlen, xlen), ierr)
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
            real (rp), dimension(:), intent(in) :: x, y
            integer (ip), intent(in) :: BCx, BCy

            ! outputs
            integer (ip), intent(out) :: ierr

            ! runtime
            type (coo_matrix) :: DDx_matrix, DDy_matrix


            ! END DECLARATIONS

            ! BEGIN FUNCTION
            DDx_matrix = double_diff_FD(x, BCx, ierr)
            DDy_matrix = double_diff_FD(y, BCy, ierr)

            DDx_matrix = sparse_kron(DDx_matrix, speye(size(y, 1_ip, KIND=ip), ierr), ierr)
            DDy_matrix = sparse_kron(speye(size(x, 1_ip, KIND=ip), ierr), DDy_matrix, ierr)
            call sparse_add(DDx_matrix, DDy_matrix, ierr)
            laplacian2D = DDx_matrix
        end function laplacian2D

        type (coo_matrix) function laplacian3D(x, y, z, BCx, BCy, BCz, ierr)
            !
            ! generates a sparse laplacian matrix
            !
            ! inputs:
            !   x, y, z: vector of equispaced grid points
            !   BCx, BCy, BCz: boundary conditions for x and y directions (1 if periodic)
            !
            ! outputs:
            !   ierr: error flag
            !   return value: laplacian matrix in sparse coo form
            !
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            real (rp), dimension(:), intent(in) :: x, y, z
            integer (ip), intent(in) :: BCx, BCy, BCz

            ! outputs
            integer (ip), intent(out) :: ierr

            ! runtime
            type (coo_matrix) :: DDx_matrix, DDy_matrix, DDz_matrix, temp


            ! END DECLARATIONS

            ! BEGIN FUNCTION
            DDx_matrix = double_diff_FD(x, BCx, ierr)
            DDx_matrix = sparse_kron(DDx_matrix, speye(size(y, 1_ip, KIND=ip), ierr), ierr)
            DDx_matrix = sparse_kron(speye(size(z, 1_ip, KIND=ip), ierr), DDx_matrix, ierr)


            DDy_matrix = double_diff_FD(y, BCy, ierr)
            temp = sparse_kron(speye(size(z, 1_ip, KIND=ip), ierr), speye(size(x, 1_ip, KIND=ip), ierr), ierr)
            DDy_matrix = sparse_kron(temp, DDy_matrix, ierr)

            DDz_matrix = double_diff_FD(z, BCz, ierr)
            temp = sparse_kron(speye(size(x, 1_ip, KIND=ip), ierr), speye(size(y, 1_ip, KIND=ip), ierr), ierr)
            DDz_matrix = sparse_kron(DDz_matrix, temp, ierr)

            call sparse_add(DDx_matrix, DDy_matrix, ierr)
            call sparse_add(DDx_matrix, DDz_matrix, ierr)
            laplacian3D = DDx_matrix
        end function laplacian3D

end module diff
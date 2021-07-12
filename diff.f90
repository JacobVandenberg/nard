! Created by  on 6/7/21.

module diff
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

        subroutine double_diff_FD(x, DDx_vals, DDx_indx, DDx_jndx, ierr)
            ! returns a sparse double derivative matrix
            ! input:
            !   N: number of gridpoints
            !   domsize: size of domain
            ! output:
            !   DDx_vals, DDx_indx, DDx_jndx: sparse double derivative matrix
            !
            implicit none
            real (kind = 8) :: dx
            real (kind = 8), dimension(:), intent(in) :: x
            real (kind = 8), intent(out), dimension(:), allocatable :: DDx_vals
            integer, intent(out), dimension(:), allocatable :: DDx_indx, DDx_jndx
            integer, intent(out) :: ierr
            integer :: xlen

            if (size(x, 1) < 2) then
                ierr = -1
                return
            end if

            dx = x(2) - x(1)

            call sparse_tridiag(dble(-2)/ (dx*dx), dble(1)/(dx*dx), dble(1)/(dx*dx), size(x, 1), DDx_vals, DDx_indx, DDx_jndx, ierr)
            return
        end subroutine double_diff_FD

end module diff
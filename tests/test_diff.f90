! Created by  on 8/7/21.

program test_diff
    use diff
    use ogpf
    implicit none
    call test_tridiag
    call test_sparse_tridiag
    Print *, 'test_double_diff_FD'
    call test_double_diff_FD
    Print *, 'test_laplacian2D'
    call test_laplacian2D
    contains
        subroutine test_tridiag
            implicit none
            integer :: N = 4, ierr
            real (kind = 8), dimension(:, :), allocatable :: result

            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate(result)

            N = 0
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate(result)

            N = 1
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate (result)

            N = 2
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate (result)
        end subroutine test_tridiag

        subroutine test_sparse_tridiag
            implicit none
            integer :: N, ierr
            real(kind=8), allocatable, dimension(:) :: val
            integer, dimension(:), allocatable :: indx, jndx

            N = 5
            call sparse_tridiag(dble(-2), dble(1), dble(-1), N, val, indx, jndx, ierr)

            Print *, val
            Print *, indx
            Print *, jndx

            deallocate(val, indx, jndx)


        end subroutine test_sparse_tridiag

        subroutine test_double_diff_FD
            implicit none
            type (coo_matrix) :: mat
            real (kind=8), dimension(:), allocatable :: x
            integer :: ierr

            x = linspace(dble(0), dble(4), 5)
            mat = double_diff_FD(x, 1, ierr)
            Print *, mat%vals
            deallocate(x)

            x = linspace(dble(0), dble(4), 5)
            mat = double_diff_FD(x, 0, ierr)
            Print *, mat%vals

        end subroutine test_double_diff_FD

        subroutine test_laplacian2D
            implicit none
            type (coo_matrix) :: mat
            type (gpf) :: gp
            real (kind=8), dimension(:), allocatable :: x, y
            real (kind=8), dimension(:, :), allocatable :: xx, yy, Luu, uu_plot
            integer :: ierr, N
            N = 50
            x = linspace(dble(-1), dble(1), N)
            y = linspace(dble(-1), dble(1), N)

            call meshgrid(xx, yy, x, y, ierr)

            mat = laplacian2D(x, y, 0, 0, ierr)
            allocate(Luu(N * N, 1), uu_plot(N, N))

            Print *, mat%n, mat%m

            call coo_multiply(mat, reshape(cos(3.14159 * xx) + cos(3.14159 * yy), (/N*N, 1/)), Luu)
            uu_plot = reshape(Luu, (/N, N/))
            call gp%surf(xx, yy, uu_plot)

        end subroutine test_laplacian2D
end program test_diff
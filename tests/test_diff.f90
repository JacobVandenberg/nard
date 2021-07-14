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
            integer (ip) :: N = 4_ip, ierr
            real (rp), dimension(:, :), allocatable :: result

            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0_ip) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate(result)

            N = 0_ip
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0_ip) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate(result)

            N = 1_ip
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0_ip) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate (result)

            N = 2_ip
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0_ip) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate (result)
        end subroutine test_tridiag

        subroutine test_sparse_tridiag
            implicit none
            integer (ip) :: N, ierr
            real(rp), allocatable, dimension(:) :: val
            integer (ip), dimension(:), allocatable :: indx, jndx

            N = 5_ip
            call sparse_tridiag(dble(-2), dble(1), dble(-1), N, val, indx, jndx, ierr)

            Print *, val
            Print *, indx
            Print *, jndx

            deallocate(val, indx, jndx)


        end subroutine test_sparse_tridiag

        subroutine test_double_diff_FD
            implicit none
            type (coo_matrix) :: mat
            real (rp), dimension(:), allocatable :: x
            integer (ip) :: ierr

            x = linspace(dble(0), dble(4), 5_ip)
            mat = double_diff_FD(x, 1_ip, ierr)
            Print *, mat%vals
            deallocate(x)

            x = linspace(dble(0), dble(4), 5_ip)
            mat = double_diff_FD(x, 0_ip, ierr)
            Print *, mat%vals

        end subroutine test_double_diff_FD

        subroutine test_laplacian2D
            implicit none
            type (coo_matrix) :: mat
            type (gpf) :: gp
            real (rp), dimension(:), allocatable :: x, y
            real (rp), dimension(:, :), allocatable :: xx, yy, Luu, uu_plot
            integer (ip) :: ierr, N
            N = 50_ip
            x = linspace(dble(-1), dble(1), N)
            y = linspace(dble(-1), dble(1), N)

            call meshgrid(xx, yy, x, y, ierr)

            mat = laplacian2D(x, y, 0_ip, 0_ip, ierr)
            allocate(Luu(N * N, 1), uu_plot(N, N))

            Print *, mat%n, mat%m

            call coo_multiply(mat, reshape(cos(3.14159 * xx) + cos(3.14159 * yy), (/N*N, 1_ip/)), Luu)
            uu_plot = reshape(Luu, (/N, N/))
            call gp%surf(xx, yy, uu_plot)

        end subroutine test_laplacian2D
end program test_diff